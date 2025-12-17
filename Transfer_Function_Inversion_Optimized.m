
%% 3B) Reconstruct Road Profile in the Frequency Domain (Transfer-Function Inversion)
% -------------------------------------------------------------------------
% PURPOSE:
% Reconstructs the road profile by dividing the sprung acceleration spectrum 
% by the analytical quarter-car transfer function, using regularization to 
% stabilize inversion at anti-resonances.
%
% CORE LOGIC:
% 1. Load: Reads `recon_cfg.mat` (accel a_s, speed V, vehicle params).
% 2. FFT: Computes Fourier transform of acceleration: A(f) = FFT{a_s(t)}.
% 3. Transfer Function: Builds analytical H(s) = Accel(s) / Road(s).
% 4. Inversion: Solves R(f) = A(f) ./ H(f) using Tikhonov or Threshold regularization.
% 5. IFFT & Map: Transforms back to time domain z_r(t), then resamples to spatial z_r(x).
% 6. Output: Saves standard spatial struct `sp` to '03B_TF Inversion/recon_spatial.mat'.
% -------------------------------------------------------------------------

close all; clc;
fprintf('\n=== PATH B: Direct TF inversion ===\n');

%% ---------------------------- HYPERPARAMETERS -----------------------------
% --- Reconstruction / inversion ---
use_tikhonov          = true;      
% Adaptive Regularization Settings
use_adaptive_lambda   = true;                  % Auto-tune lambda using L-Curve?
tikhonov_lambda_init  = 1e-5;                  % Initial guess / fallback
lambda_scan_range     = logspace(-8, -1, 50);  % Range to scan for optimal lambda

% Signal Processing Improvements
tukey_taper_ratio     = 0.05;      % Taper 5% of edges to reduce spectral leakage
f_min_enforce         = 0.05;       % Post-inversion High-Pass cutoff [Hz] (removes drift)

drop_dB               = 40;        % keep freqs within ~drop_dB dB of peak |H|
f_min                 = 0.1;       % min inversion frequency [Hz]
f_max_cap             = 90;        % max inversion frequency cap [Hz]; actual f_max = min(f_max_cap, 0.95*fs/2)

% IRI-relevant band & lambda shaping
f_iri_min             = 0.1;       % lower edge of IRI-relevant band [Hz]
f_iri_max             = 50;        % upper edge of IRI-relevant band [Hz]
lambda_iri_factor     = 1e-2;      % lambda multiplier inside IRI band (milder reg.)
lambda_outband_factor = 1e+1;      % lambda multiplier outside IRI band (stronger reg.)

% --- Spatial mapping ---
spatial_grid_policy   = 'ground_truth'; % {'ground_truth','uniform'}

% --- Output / publishing ---
publish_to_shared     = true;
out_root_dir          = '00_Outputs';
out_dir_pathB         = fullfile(out_root_dir, '03_ReconstructionCore', 'PathB_FrequencyDomain');

% --- Plot controls ---
make_plots            = true;
save_plots            = true;      % write PNG + FIG under out_dir_pathB/figs
fig_format            = 'png';     % {'png','tiff','pdf','png'} (pdf not for desktop transparency)
fig_dpi               = 300;       % export resolution
use_latex_labels      = false;     % set true if you want LaTeX (requires proper setup)

% Plot style
baseFont              = 11;
axesLW                = 0.9;       % (used as a reference; actual value set in applyTheme)
lineLW                = 1.6;
useMinorGrid          = true;

%% ----------------------------- OUTPUT FOLDERS -----------------------------
fig_dir = fullfile(out_dir_pathB,'figs');
if ~exist(out_dir_pathB,'dir')
    mkdir(out_dir_pathB);
end
if make_plots && save_plots && ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% -------------------------- INPUT AUTO-DISCOVERY --------------------------
cfg_candidates = { ...
    fullfile(out_root_dir, '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};

cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file')
        cfg_path = cfg_candidates{k};
        break;
    end
end
assert(~isempty(cfg_path), 'PATH B: recon_cfg.mat not found. Run prepare_reconstruction_inputs.m first.');

S = load(cfg_path);  assert(isfield(S,'cfg'), 'PATH B: recon_cfg.mat missing cfg struct.');
cfg = S.cfg;

t = cfg.t(:); fs = cfg.fs; V = cfg.V;

accel_field_order = {'az_s','a_z_s','a_s','azs','a_sprung'};
az_s = [];
for f = accel_field_order
    if isfield(cfg, f{1})
        az_s = cfg.(f{1})(:);
        break;
    end
end
assert(~isempty(az_s),'PATH B: Could not locate sprung acceleration in cfg (tried %s).', strjoin(accel_field_order,', '));

ms = cfg.ms;
mu = cfg.mu;
ks = cfg.ks;
cs = cfg.cs;
kt = cfg.kt;
ct = 0;

if isfield(cfg,'ct') && ~isempty(cfg.ct)
    ct = cfg.ct;
end

if isfield(cfg,'x_t') && ~isempty(cfg.x_t)
    x_t = cfg.x_t(:);
else
    x_t = (t - t(1)) * V;
end

N = numel(t);
fprintf('Inputs: N=%d, fs=%.3f Hz, V=%.3f m/s\n', N, fs, V);

%% -------------------- 1) FFT of sprung acceleration ----------------------
% Remove DC bias first
az_s = az_s - mean(az_s);

% Apply Tukey Window to reduce spectral leakage at edges
win_taper = tukeywin(N, tukey_taper_ratio);
az_s_win  = az_s .* win_taper;

% FFT (Optional: Zero-pad to next power of 2 * 2 for smoother spectrum)
N_fft   = 2 * 2^nextpow2(N); 
Accel_f = fft(az_s_win, N_fft);
freq    = (0:N_fft-1)' * fs / N_fft;
s       = 1i * 2*pi*freq;

%% --------------------- 2) Quarter-car transfer H(s) ----------------------
% H(s): A_s(s)/Z_r(s) — from road displacement to sprung-body acceleration
num_H = (s.^2) .* (kt*cs.*s + kt*ks);
den_H = (ms*mu).*s.^4 ...
      + (ms*ct + ms*cs + mu*cs).*s.^3 ...
      + (ms*kt + ms*ks + mu*ks + cs*ct).*s.^2 ...
      + (ks*ct + kt*cs).*s ...
      + (ks*kt);
Hf = num_H ./ den_H;

% ---- Adaptive hard threshold: X dB below peak |H| ----
Hmag_full  = abs(Hf);
fnyq       = fs/2;
band       = (freq > 0) & (freq < 0.95*fnyq);   % ignore DC and extreme top
if ~any(band)
    band = freq > 0;                            % fallback: just ignore DC
end
Hmax = max(Hmag_full(band));
if Hmax <= 0 || ~isfinite(Hmax)
    inversion_threshold = 0;                    % fallback: effectively no hard cut
else
    inversion_threshold = Hmax * 10^(-drop_dB/20);
end

% Natural frequencies (from polynomial denominator coefficients)
poly_den = [ms*mu, (ms*ct + ms*cs + mu*cs), ...
            (ms*kt + ms*ks + mu*ks + cs*ct), ...
            (ks*ct + kt*cs), ks*kt];
poles = roots(poly_den);                               % complex poles (rad/s)
wn_hz = sort(unique(round(abs(imag(poles))/2/pi, 4))); % indicative resonances (Hz)

%% --------------------- 3) Regularized frequency inversion ----------------
f_max = min(f_max_cap, 0.95*(fs/2));
band_mask = (freq >= f_min) & (freq <= f_max);

if use_tikhonov
    fprintf('Inversion: Tikhonov regularization on [%.2f, %.2f] Hz\n', f_min, f_max);
    
    Hmag2 = abs(Hf).^2;
    
    % --- ADAPTIVE L-CURVE OPTIMIZATION ---
    if use_adaptive_lambda
        fprintf('  -> Running L-Curve optimization... ');
        
        % Extract valid band data for optimization (faster)
        Hf_b = Hf(band_mask);
        Af_b = Accel_f(band_mask);
        H2_b = Hmag2(band_mask);
        
        eta = zeros(size(lambda_scan_range)); % Solution Norm ||R||
        rho = zeros(size(lambda_scan_range)); % Residual Norm ||HR - A||
        
        for k = 1:length(lambda_scan_range)
            lam = lambda_scan_range(k);
            % Tikhonov estimate on band
            R_est_b = (conj(Hf_b) ./ (H2_b + lam)) .* Af_b;
            
            % L-Curve norms
            eta(k) = norm(R_est_b);                  
            rho(k) = norm((Hf_b .* R_est_b) - Af_b); 
        end
        
        % Find L-Curve "Corner" (closest point to origin in log-log space)
        log_eta = log10(eta); 
        log_rho = log10(rho);
        % Normalize to [0,1]
        n_eta = (log_eta - min(log_eta)) / (max(log_eta) - min(log_eta));
        n_rho = (log_rho - min(log_rho)) / (max(log_rho) - min(log_rho));
        
        [~, best_idx] = min(n_eta.^2 + n_rho.^2);
        tikhonov_lambda = lambda_scan_range(best_idx);
        
        fprintf('Optimal Lambda = %.2e\n', tikhonov_lambda);
        
        % Optional: Plot L-Curve for debugging
        if make_plots
            figure('Name','L-Curve','Visible','off'); loglog(rho, eta, 'b-'); hold on;
            loglog(rho(best_idx), eta(best_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
            title(sprintf('L-Curve (Opt \\lambda = %.2e)', tikhonov_lambda));
            xlabel('Residual Norm ||HR - A||'); ylabel('Solution Norm ||R||'); grid on;
            exportgraphics(gcf, fullfile(fig_dir, 'B00_L_Curve.png'));
        end
    else
        tikhonov_lambda = tikhonov_lambda_init;
    end

    % --- APPLY INVERSION ---
    Road_f = zeros(size(Accel_f));
    
    % Apply IRI-band weighting (optional preservation of your original logic)
    lambda_vec = tikhonov_lambda * ones(size(freq));
    % If you want to keep the IRI specific scaling:
    if exist('f_iri_min','var')
        iri_band = (freq >= max(f_iri_min, f_min)) & (freq <= min(f_iri_max, f_max));
        lambda_vec(iri_band) = tikhonov_lambda * lambda_iri_factor;
        lambda_vec(band_mask & ~iri_band) = tikhonov_lambda * lambda_outband_factor;
    end

    reg = Hmag2 + lambda_vec;
    Road_f(band_mask) = (conj(Hf(band_mask)) ./ reg(band_mask)) .* Accel_f(band_mask);
    
    valid = band_mask;
else
    % (Keep your existing Hard Threshold logic here if you want)
end

%% --------------------- 4) Back to time, then to spatial ------------------
zr_t_padded = ifft(Road_f, 'symmetric');
zr_t = zr_t_padded(1:N);

% --- CRITICAL IMPROVEMENT: Post-Inversion Detrending ---
% Frequency domain inversion leaves low-freq drift. We must remove it.
% 1. Remove Linear Trend
zr_t = detrend(zr_t, 'linear');

% 2. High-Pass Filter (remove integration drift < f_min)
if f_min_enforce > 0
    [b_hp, a_hp] = butter(2, f_min_enforce / (fs/2), 'high');
    zr_t = filtfilt(b_hp, a_hp, zr_t);
end

switch lower(spatial_grid_policy)
    case 'ground_truth'
        if isfield(cfg,'ground') && isfield(cfg.ground,'x')
            x_grid = cfg.ground.x(:);
        else
            warning('PATH B: ground truth grid missing; falling back to UNIFORM.');
            x_grid = (t - t(1))*V;
        end
    case 'uniform'
        x_grid = (t - t(1))*V;
    otherwise
        error('PATH B: Unknown spatial_grid_policy: %s', spatial_grid_policy);
end

zr_x_raw  = interp1(x_t, zr_t, x_grid, 'linear', 'extrap');
zr_x_raw(~isfinite(zr_x_raw)) = 0;
zr_x_filt = zr_x_raw;

dx         = mean(diff(x_grid),'omitnan');
fs_spatial = 1 / max(dx, eps);

sp = struct();
sp.x              = x_grid;
sp.zr_x_raw       = zr_x_raw;
sp.zr_x_filt      = zr_x_filt;
sp.dx             = dx;
sp.fs_spatial     = fs_spatial;
sp.nc_spatial     = NaN;
sp.lp_order       = NaN;
sp.filter_applied = false;
sp.meta = struct( ...
    'recon_method',         'Path B: Direct TF Inversion', ...
    'inversion_threshold',   inversion_threshold, ...
    'use_tikhonov',          use_tikhonov, ...
    'tikhonov_lambda',       tikhonov_lambda, ...
    'grid_policy',           spatial_grid_policy, ...
    'fs',                    fs, ...
    'V',                     V);

% Save Path B artifact
out_mat_pathB = fullfile(out_dir_pathB,'recon_spatial.mat');
save(out_mat_pathB, 'sp','-v7.3');
fprintf('PATH B: saved %s\n', out_mat_pathB);

% Optional publish for downstream compatibility
if publish_to_shared
    share_dir = fullfile(out_root_dir, '03_ReconstructionCore', 'SpatialMapping');
    if ~exist(share_dir,'dir')
        mkdir(share_dir);
    end
    copyfile(out_mat_pathB, fullfile(share_dir,'recon_spatial.mat'));
    fprintf('PATH B: also published to %s\n', share_dir);
end

%% ------------------------------ 5) PLOTS ----------------------------------
if make_plots
    % ---------- Global style ----------
    if use_latex_labels
        set(0,'defaultTextInterpreter','latex');
        set(0,'defaultLegendInterpreter','latex');
        set(0,'defaultAxesTickLabelInterpreter','latex');
    end

    set(0,'defaultAxesFontName','Helvetica',...
          'defaultAxesFontSize',baseFont,...
          'defaultLineLineWidth',lineLW,...
          'DefaultAxesBox','off');   % cleaner, modern look

    % Color palette (colorblind-friendly-ish, slightly softened)
    C.black = [0.10 0.10 0.10];
    C.blue  = [0.00 0.45 0.74];
    C.red   = [0.80 0.25 0.25];
    C.gray  = [0.55 0.55 0.60];
    C.green = [0.20 0.62 0.20];
    C.purp  = [0.42 0.24 0.60];
    C.iri   = [0.40 0.75 0.65];  % highlight IRI band

    % Frequency vectors excluding DC (for log x)
    k1   = max(2, find(freq>0,1,'first'));  % avoid f = 0
    fplt = freq(k1:end);
    Hmag = abs(Hf(k1:end));
    Hdb  = mag2db(max(Hmag, eps));
    Hphi = unwrap(angle(Hf(k1:end)))*180/pi;
    fnyq = fs/2;

    % Helper for IRI band shading (if those vars exist)
    has_iri_band = exist('f_iri_min','var') && exist('f_iri_max','var');

    % ======================================================================
    % Fig 1: Bode (|H| in dB and phase)
    % ======================================================================
    fig1 = figure('Name','Path B — Transfer Function',...
                  'Color','w','Units','pixels');
    fig1.Position(3:4) = [980 640];
    tlo1 = tiledlayout(fig1,2,1,...
        'TileSpacing','compact','Padding','compact');

    % Magnitude
    ax1 = nexttile(tlo1,1); hold(ax1,'on');
    set(ax1,'XScale','log');

    % Magnitude curve
    hMag = plot(ax1, fplt, Hdb, ...
        'Color', C.black, 'LineWidth', 1.8, ...
        'DisplayName','|H(f)|');

    % Optional shaded IRI band
    if use_tikhonov && has_iri_band
        iri_min = max(f_iri_min, fplt(1));
        iri_max = min(f_iri_max, fnyq);
        if iri_max > iri_min
            yl = ylim(ax1);
            patch(ax1, ...
                [iri_min iri_max iri_max iri_min], ...
                [yl(1)   yl(1)   yl(2)   yl(2)], ...
                C.iri, 'FaceAlpha',0.06, ...
                'EdgeColor','none', ...
                'HandleVisibility','off');
            uistack(hMag,'top');
        end
    end

    % Threshold line
    thrLine = yline(ax1, mag2db(inversion_threshold), '--', ...
        'Color', C.red, 'LineWidth',1.1, ...
        'DisplayName','|H| threshold'); %#ok<NASGU>

    % Masked regions for hard threshold
    if ~use_tikhonov
        drawnow;   % ensure ylim exists
        yl1 = ylim(ax1);
        shade_invalid_regions(ax1, fplt, valid(k1:end), yl1, C.red, 0.08);
    end

    % Nyquist line
    xline(ax1, fnyq, ':', 'Nyquist', ...
        'Color', C.gray, 'LineWidth',1.0, ...
        'LabelVerticalAlignment','top', ...
        'LabelHorizontalAlignment','left', ...
        'HandleVisibility','off');

    % Resonance markers
    for k = 1:numel(wn_hz)
        if wn_hz(k) > fplt(1) && wn_hz(k) < fnyq
            xline(ax1, wn_hz(k), ':', ...
                sprintf('%.2f Hz', wn_hz(k)), ...
                'Color', C.gray, 'LineWidth',0.8, ...
                'LabelVerticalAlignment','bottom', ...
                'LabelHorizontalAlignment','left', ...
                'HandleVisibility','off');
        end
    end

    applyTheme(ax1, useMinorGrid);
    ylabel(ax1,'Magnitude |H(f)| [dB]');
    title(ax1,'Vehicle transfer function H(f): road \rightarrow sprung acceleration');

    % Legend slightly inset and compact
    niceLegend(ax1,'southwest');

    % Phase
    ax2 = nexttile(tlo1,2); hold(ax2,'on');
    set(ax2,'XScale','log');
    plot(ax2, fplt, Hphi, ...
        'Color', C.blue, 'LineWidth',1.6, ...
        'DisplayName','Phase');

    xline(ax2, fnyq, ':', ...
        'Color', C.gray, 'LineWidth',1.0, ...
        'HandleVisibility','off');

    for k = 1:numel(wn_hz)
        if wn_hz(k) > fplt(1) && wn_hz(k) < fnyq
            xline(ax2, wn_hz(k), ':', ...
                'Color', C.gray, 'LineWidth',0.8, ...
                'HandleVisibility','off');
        end
    end

    applyTheme(ax2, useMinorGrid);
    xlabel(ax2,'Frequency [Hz]');
    ylabel(ax2,'Phase [deg]');
    title(ax2,'Phase response');
    linkaxes([ax1,ax2],'x');
    xlim(ax1,[fplt(1) fnyq*1.05]);

    % Tiny annotation: fraction of usable spectrum
    fracValid = 100*mean(valid(k1:end));
    annotation(fig1,'textbox',[.67 .91 .30 .06], ...
        'String',sprintf('Valid frequencies used: %.1f%%', fracValid), ...
        'FitBoxToText','on','EdgeColor','none','Color',C.purp, ...
        'HorizontalAlignment','right','FontSize',baseFont-1);

    export_pub(fig1, fig_dir, 'B01_TF_Bode', fig_format, fig_dpi, save_plots);

    % ======================================================================
    % Fig 2: Spectra — |A(f)| and |R(f)| (dB)
    % ======================================================================
    Af_mag = abs(Accel_f(k1:end));  Af_db = mag2db(max(Af_mag, eps));
    Rf_mag = abs(Road_f(k1:end));   Rf_db = mag2db(max(Rf_mag, eps));

    fig2 = figure('Name','Path B — Spectra',...
                  'Color','w','Units','pixels');
    fig2.Position(3:4) = [980 640];
    tlo2 = tiledlayout(fig2,2,1,...
        'TileSpacing','compact','Padding','compact');

    % Acceleration spectrum
    ax3 = nexttile(tlo2,1); hold(ax3,'on');
    set(ax3,'XScale','log');

    plot(ax3, fplt, Af_db, ...
        'Color', C.black, 'LineWidth',1.8, ...
        'DisplayName','|A(f)|  (sprung accel.)');

    xline(ax3, fnyq, ':', ...
        'Color', C.gray, 'LineWidth',1.0, ...
        'HandleVisibility','off');

    applyTheme(ax3, useMinorGrid);
    ylabel(ax3,'Amplitude [dB re 1 m/s^2]');
    title(ax3,'Measured acceleration spectrum');
    niceLegend(ax3,'southwest');

    % Reconstructed road spectrum
    ax4 = nexttile(tlo2,2); hold(ax4,'on');
    set(ax4,'XScale','log');

    plot(ax4, fplt, Rf_db, ...
        'Color', C.blue, 'LineWidth',1.8, ...
        'DisplayName','|R(f)|  (reconstructed road)');

    xline(ax4, fnyq, ':', ...
        'Color', C.gray, 'LineWidth',1.0, ...
        'HandleVisibility','off');

    if ~use_tikhonov
        drawnow;
        yl4 = ylim(ax4);
        shade_invalid_regions(ax4, fplt, valid(k1:end), yl4, C.red, 0.08);
    end

    applyTheme(ax4, useMinorGrid);
    xlabel(ax4,'Frequency [Hz]');
    ylabel(ax4,'Amplitude [dB re 1 m]');
    title(ax4,'Reconstructed road spectrum');
    linkaxes([ax3,ax4],'x');
    xlim(ax3,[fplt(1) fnyq*1.05]);
    niceLegend(ax4,'southwest');

    export_pub(fig2, fig_dir, 'B02_Spectra', fig_format, fig_dpi, save_plots);

    % ======================================================================
    % Fig 3: Spatial profile vs ground + error panel
    % ======================================================================
    if isfield(cfg,'ground') && isfield(cfg.ground,'x') && isfield(cfg.ground,'zr')
        xg = cfg.ground.x(:);
        zg = cfg.ground.zr(:);
        [x_common, zr_est_resampled, zg_resampled] = ...
            align_on_common_grid(sp.x, sp.zr_x_filt, xg, zg);
        err = zr_est_resampled - zg_resampled;

        % Metrics
        rmse = sqrt(mean(err.^2,'omitnan'));
        mae  = mean(abs(err),'omitnan');

        vp = isfinite(zr_est_resampled) & isfinite(zg_resampled);
        if nnz(vp) >= 2 && std(zr_est_resampled(vp)) > 0 && std(zg_resampled(vp)) > 0
            rho = corr(zr_est_resampled(vp), zg_resampled(vp));
        else
            rho = NaN;
        end

        % Symmetric error limits (robust)
        emax = robust_max(abs(err));
        emax = max(emax, eps);
        elims = 1.05*[-emax, emax];

        fig3 = figure('Name','Path B — Spatial Comparison',...
                      'Color','w','Units','pixels');
        fig3.Position(3:4) = [980 640];
        tlo3 = tiledlayout(fig3,2,1,...
            'TileSpacing','compact','Padding','compact');

        % Top: ground vs recon
        ax5 = nexttile(tlo3,1); hold(ax5,'on');
        plot(ax5, xg, zg, '-', ...
            'Color', C.gray, 'LineWidth',1.5, ...
            'DisplayName','Ground truth');
        plot(ax5, sp.x, sp.zr_x_filt, '-', ...
            'Color', C.blue, 'LineWidth',1.8, ...
            'DisplayName','TF inversion');

        applyTheme(ax5, useMinorGrid);
        ylabel(ax5,'Elevation z_r(x) [m]');
        title(ax5, sprintf(['Spatial profile  (RMSE = %.3g m, ', ...
                            'MAE = %.3g m, \\rho = %.3f)'], ...
                            rmse, mae, rho));
        niceLegend(ax5,'best');

        % Bottom: error
        ax6 = nexttile(tlo3,2); hold(ax6,'on');
        plot(ax6, x_common, err, '-', ...
            'Color', C.red, 'LineWidth',1.4, ...
            'DisplayName','Error (recon − truth)');

        yline(ax6, 0, ':', ...
            'Color', C.gray, 'LineWidth',1.0, ...
            'HandleVisibility','off');

        applyTheme(ax6, useMinorGrid);
        xlabel(ax6,'Distance x [m]');
        ylabel(ax6,'Error [m]');
        % ylim(ax6, elims);
        title(ax6,'Reconstruction error (symmetric limits)');
        linkaxes([ax5,ax6],'x');
        niceLegend(ax6,'best');

        export_pub(fig3, fig_dir, 'B03_Spatial_Comparison', fig_format, fig_dpi, save_plots);
    end

    % ======================================================================
    % Fig 4: Time-domain sanity checks
    % ======================================================================
    fig4 = figure('Name','Path B — Time Signals',...
                  'Color','w','Units','pixels');
    fig4.Position(3:4) = [980 640];
    tlo4 = tiledlayout(fig4,2,1,...
        'TileSpacing','compact','Padding','compact');

    % Sprung acceleration
    ax7 = nexttile(tlo4,1); hold(ax7,'on');
    plot(ax7, t, az_s, ...
        'Color', C.black, 'LineWidth',1.5, ...
        'DisplayName','a_s(t)');
    applyTheme(ax7, useMinorGrid);
    xlim(ax7, [t(1) t(end)]);
    ylabel(ax7,'Sprung accel. a_s(t) [m/s^2]');
    title(ax7,'Measured sprung acceleration (time domain)');
    niceLegend(ax7,'best');

    % Reconstructed road vs distance
    ax8 = nexttile(tlo4,2); hold(ax8,'on');
    x_dist = (t - t(1))*V;
    plot(ax8, x_dist, zr_t, ...
        'Color', C.green, 'LineWidth',1.6, ...
        'DisplayName','z_r(t) mapped to x');
    applyTheme(ax8, useMinorGrid);
    xlim(ax8, [min(x_dist) max(x_dist)]);
    xlabel(ax8,'Distance x = V(t - t_0) [m]');
    ylabel(ax8,'Reconstructed z_r [m]');
    title(ax8,'Reconstructed road before spatial resampling');
    niceLegend(ax8,'best');

    export_pub(fig4, fig_dir, 'B04_TimeSignals', fig_format, fig_dpi, save_plots);
end

fprintf('=== PATH B: done ===\n\n');

%% ------------------------------ HELPERS -----------------------------------
function shade_invalid_regions(ax, f, validMask, ylims, color, alphaVal)
% Shades regions where validMask==false along the frequency axis
    if all(validMask)
        return;
    end
    invalid = ~validMask(:);
    % Find contiguous invalid segments
    d = diff([false; invalid; false]);
    starts = find(d==1);  ends = find(d==-1)-1;
    for i = 1:numel(starts)
        fx = [f(starts(i)) f(ends(i)) f(ends(i)) f(starts(i))];
        fy = [ylims(1) ylims(1) ylims(2) ylims(2)];
        patch('Parent',ax,'XData',fx,'YData',fy,'FaceColor',color,'FaceAlpha',alphaVal, ...
              'EdgeColor','none','HandleVisibility','off');
    end
end

function maybe_export(fig, outdir, basename, fmt, dpi, doSave)
% Exports figure as PNG/TIFF/PDF + MATLAB FIG if requested
    if ~doSave, return; end
    fp = fullfile(outdir, basename);
    switch lower(fmt)
        case 'png'
            print(fig, [fp '.png'], '-dpng', sprintf('-r%d',dpi));
        case 'tiff'
            print(fig, [fp '.tif'], '-dtiff', sprintf('-r%d',dpi));
        case 'pdf'
            set(fig,'PaperPositionMode','auto'); print(fig, [fp '.pdf'], '-dpdf','-vector');
        otherwise
            print(fig, [fp '.png'], '-dpng', sprintf('-r%d',dpi));
    end
    try savefig(fig, [fp '.fig']);
    catch
    end
end

function [x_common, y1i, y2i] = align_on_common_grid(x1,y1,x2,y2)
% Safe linear alignment on a common grid (for error plots/metrics)
    xmin = max(min(x1), min(x2));
    xmax = min(max(x1), max(x2));
    if xmax <= xmin
        x_common = [];
        y1i = [];
        y2i = [];
        return;
    end
    % Choose a reasonable step from both inputs
    dx = median([median(diff(x1),'omitnan'), median(diff(x2),'omitnan')],'omitnan');
    if ~isfinite(dx) || dx <= 0
        dx = (xmax - xmin) / max(1000, numel(x1) + numel(x2));
    end
    x_common = (xmin:dx:xmax).';
    y1i = interp1(x1,y1,x_common,'linear','extrap');
    y2i = interp1(x2,y2,x_common,'linear','extrap');
end

function out = robust_max(x)
% Returns a robust upper scale (98th percentile) to avoid wild y-lims
    x = x(isfinite(x));
    if isempty(x)
        out = 0;
        return;
    end
    out = prctile(x,98);
end

function export_pub(fig, outdir, basename, fmt, dpi, doSave)
% exportgraphics when available (vector PDF), fallback to print
    if ~doSave
        return;
    end
    fp = fullfile(outdir, basename);
    try
        switch lower(fmt)
            case 'png'
                % No transparency request -> no warning
                exportgraphics(fig, [fp '.png'], 'Resolution', dpi);
            case 'tiff'
                exportgraphics(fig, [fp '.tif'], 'Resolution', dpi);
            case 'pdf'
                exportgraphics(fig, [fp '.pdf'], 'ContentType','vector');
            otherwise
                exportgraphics(fig, [fp '.png'], 'Resolution', dpi);
        end
    catch
        % Fallback for older MATLAB
        switch lower(fmt)
            case 'png'
                print(fig, [fp '.png'], '-dpng', sprintf('-r%d',dpi));
            case 'tiff'
                print(fig, [fp '.tif'], '-dtiff', sprintf('-r%d',dpi));
            case 'pdf'
                set(fig,'PaperPositionMode','auto'); print(fig, [fp '.pdf'], '-dpdf','-vector');
            otherwise
                print(fig, [fp '.png'], '-dpng', sprintf('-r%d',dpi));
        end
    end
    try savefig(fig, [fp '.fig']);
    catch
    end
end

function applyTheme(ax, useMinor)
% Consistent axis styling for all Path B figures.
    if nargin < 2
        useMinor = true;
    end
    set(ax, ...
        'Box','on', ...
        'LineWidth',0.9, ...
        'FontName','Helvetica', ...
        'FontSize',11, ...
        'XGrid','on','YGrid','on');

    if useMinor
        set(ax,'XMinorGrid','on','YMinorGrid','on');
    else
        set(ax,'XMinorGrid','off','YMinorGrid','off');
    end
end

function niceLegend(ax, loc)
% Compact legend with consistent style.
    leg = legend(ax,'Location',loc);
    if ~isempty(leg) && isvalid(leg)
        set(leg, 'Box','off', 'ItemTokenSize',[10 8], 'FontSize',9);
    end
end