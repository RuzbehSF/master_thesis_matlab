
%% 5) Calculate Reconstruction Quality and Confidence Intervals
% -------------------------------------------------------------------------
% PURPOSE:
% Assesses the quality and stability of the reconstructed road profile
% using a bootstrap methodology. It generates confidence intervals (CI) for
% the profile and calculates segment-by-segment quality metrics.
%
% CORE LOGIC (How it works):
% 1.  Loads the final spatial profile ('recon_spatial.mat') and the main
%     configuration file ('recon_cfg.mat'), which contains the original
%     clean acceleration signal.
%
% 2.  Estimates the level of sensor noise present in the system, either by
%     using a fixed value or by analyzing the high-frequency content of the
%     original acceleration signal's Power Spectral Density (PSD).
%
% 3.  Runs a bootstrap analysis. This involves a loop that runs multiple
%     times (e.g., 20-50 iterations). In each iteration, it adds a new,
%     randomly generated noise signal to the clean acceleration data.
%
% 4.  For each noisy acceleration signal, it re-runs the ENTIRE
%     reconstruction pipeline (time-domain inversion and spatial mapping)
%     to produce a new version of the final road profile.
%
% 5.  After collecting all the bootstrapped profile versions, it performs
%     statistical analysis at each point along the road to calculate the
%     mean profile, the standard deviation (a measure of stability), and
%     the 95% confidence interval.
%
% 6.  Calculates quality metrics (e.g., Signal-to-Noise Ratio, stability)
%     for fixed-length segments of the road.
%
% 7.  Saves all CI and quality metrics to a new .mat file.
%
% INPUTS:
% - '02_Reconstruction Inputs/recon_cfg.mat': For the original acceleration
%   signal and vehicle parameters.
% - '04_Spatial Mapping/recon_spatial.mat': The baseline reconstructed
%   profile.
%
% OUTPUTS:
% - '05_Quality and CI/quality_results.mat': Contains the calculated
%   confidence intervals and segment quality metrics.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ============================== HYPERPARAMETERS ==============================
% (1) Input discovery (new locations first, then legacy)
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};
sp_candidates  = { ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'SpatialMapping', 'recon_spatial.mat'), ...
};

% (2) This script’s output folder (per-script organization)
save_folder = fullfile('00_Outputs', '04_Validation', 'QualityAndCI');
fig_dir     = fullfile(save_folder,'figs');

% (3) Which reconstructed series to place CI on
%     'filtered' -> sp.zr_x_filt (recommended); 'raw' -> sp.zr_x_raw
recon_series = 'filtered';  % 'filtered' | 'raw'

% (4) Bootstrap settings
n_boot          = 20;        % number of bootstrap replicates (20–50 is plenty)
rng_seed        = 7;         % reproducibility
noise_method    = 'hf_psd';  % 'hf_psd' | 'fixed_g'
fixed_sigma_g   = 0.01;      % g-rms if noise_method='fixed_g'
hf_band_factor  = 1.5;       % estimate noise from f >= factor * f_wh

% (5) Paper-default reverse pipeline (drift control & derivatives)
hp_disp_fc_hz   = 0.05;      % high-pass cutoff for z_s (Hz) [paper]
hp_disp_order   = 2;         % IIR order
use_sgolay      = true;      % Savitzky–Golay for derivatives
sgolay_win_sec  = 0.15;      % s
sgolay_poly     = 3;         % polynomial order

% (6) Quality metrics (segment-level)
segment_length_m   = 100;     % segment size for metrics
snr_good_db        = 20;      % SNR threshold for "good"
stab_good_mm       = 2.0;     % bootstrap std (mm) threshold for "good"
resprox_good_frac  = 0.50;    % resonance power fraction threshold for "good"

% (7) Plots
make_plots = true;
export_png = true;

%% ------------------------------- Locate & load -------------------------------
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file') == 2
        cfg_path = cfg_candidates{k};
        break;
    end
end
assert(~isempty(cfg_path), 'recon_cfg.mat not found in expected locations.');

sp_path = '';
for k = 1:numel(sp_candidates)
    if exist(sp_candidates{k},'file') == 2
        sp_path = sp_candidates{k};
        break;
    end
end
assert(~isempty(sp_path), 'recon_spatial.mat not found in expected locations.');

C = load(cfg_path);
cfg = C.cfg;
SP = load(sp_path);
sp = SP.sp;

% Ensure output folders
if ~exist(save_folder,'dir')
    mkdir(save_folder);
end
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% ------------------------------ Prepare series ------------------------------
x_grid = sp.x(:);
switch lower(recon_series)
    case 'filtered'
        zr_base = sp.zr_x_filt(:);
    case 'raw'
        zr_base = sp.zr_x_raw(:);
    otherwise, error('Unknown recon_series="%s"', recon_series);
end

% Essentials from cfg
assert(all(isfield(cfg, {'t','fs','V','x_t','ms','mu','ks','cs','kt','ct'})), ...
    'cfg is missing required fields (t,fs,V,x_t,ms,mu,ks,cs,kt,ct).');

t   = cfg.t(:);
fs  = cfg.fs;
V   = cfg.V;
x_t = cfg.x_t(:);

if isfield(cfg,'az_s') && ~isempty(cfg.az_s)
    az_s = cfg.az_s(:);
else
    warning('cfg.az_s is missing. Bootstrap cannot run; writing empty CI fields.');
    sp.ci = struct('method','none','N',0,'mean',[],'std',[],'lo',[],'hi',[],'series',recon_series);
    save(sp_path,'sp','-v7.3');
    fprintf('Updated (CI skipped): %s\n', sp_path);
    return;
end

% Spatial sampling
dx = mean(diff(x_grid),'omitnan');
fs_spatial = 1/max(dx,eps);

%% ------------------------ Noise level estimation ----------------------------
g0 = 9.80665;
switch lower(noise_method)
    case 'fixed_g'
        sigma_mps2 = fixed_sigma_g * g0;
        noise_info = sprintf('fixed \\sigma = %.3f g', fixed_sigma_g);

    case 'hf_psd'
        % Estimate wheel-hop (nominal) and use PSD above hf_band_factor * f_wh
        f_wh_nom = (1/(2*pi))*sqrt(cfg.kt/cfg.mu);
        [P,f] = pwelch(az_s - mean(az_s,'omitnan'), [], [], [], fs, 'onesided');
        m_hf = f >= hf_band_factor * f_wh_nom;
        if any(m_hf)
            power_hf = trapz(f(m_hf), P(m_hf));       % approximate HF power
            sigma_mps2 = sqrt(max(power_hf, eps));    % RMS accel noise
        else
            sigma_mps2 = 0.01 * g0;                   % fallback 0.01 g
        end
        noise_info = sprintf('HF-PSD est: f >= %.1f Hz, \\sigma \\approx %.3f g', ...
            hf_band_factor*f_wh_nom, sigma_mps2/g0);

    otherwise
        error('Unknown noise_method="%s"', noise_method);
end

%% --------------------------- Bootstrap reconstructions ----------------------
rng(rng_seed);
Z = zeros(numel(x_grid), n_boot);
for b = 1:n_boot
    az_noisy = az_s + sigma_mps2*randn(size(az_s));

    Z(:,b) = recon_from_accel_once( ...
        t, az_noisy, cfg, x_t, x_grid, ...
        hp_disp_fc_hz, hp_disp_order, use_sgolay, sgolay_win_sec, sgolay_poly, ...
        recon_series, sp);  % uses sp to inherit spatial LP if 'filtered'
end

% Pointwise CI & stats
mZ  = mean(Z,2);
sZ  = std(Z,0,2);
ci  = 1.96*sZ;               % 95% normal approx
loZ = mZ - ci;
hiZ = mZ + ci;

% Create a new struct for this script's output
quality_output           = struct();
quality_output.ci        = struct();
quality_output.ci.method = sprintf('bootstrap_gaussian (%s)', noise_info);
quality_output.ci.N      = n_boot;
quality_output.ci.mean   = mZ;
quality_output.ci.std    = sZ;
quality_output.ci.lo     = loZ;
quality_output.ci.hi     = hiZ;
quality_output.ci.series = recon_series;
quality_output.ci.meta.source_sp_path = sp_path; % Keep a reference to the input file

%% ------------------------- Segment quality metrics --------------------------
edges = x_grid(1):segment_length_m:x_grid(end);
if edges(end) < x_grid(end)
    edges = [edges x_grid(end)];
end
nSeg  = numel(edges)-1;

f_wh_nom = (1/(2*pi))*sqrt(cfg.kt/cfg.mu);
[PSD_all, f_all] = pwelch(az_s - mean(az_s,'omitnan'), [], [], [], fs, 'onesided');

Q = nan(nSeg, 6);
for k = 1:nSeg
    xs = edges(k); xe = edges(k+1);

    mt = (x_t >= xs) & (x_t < xe);
    mx = (x_grid >= xs) & (x_grid < xe);

    % SNR in time (rough): var in window vs noise variance estimate
    snr_db = NaN;
    if any(mt)
        sig_var = var(az_s(mt), 1, 'omitnan');  % population var
        snr_db  = 10*log10( sig_var / max(sigma_mps2^2, eps) );
    end

    % Stability (bootstrap std) in space
    stab_std = mean(sZ(mx), 'omitnan');         % [m]
    stab_std_mm = 1e3 * stab_std;

    % Resonance proximity: fraction of PSD within ±20% of f_wh
    prox = NaN;
    if ~isempty(PSD_all)
        m_band = (f_all >= 0.8*f_wh_nom) & (f_all <= 1.2*f_wh_nom);
        frac = sum(PSD_all(m_band)) / max(sum(PSD_all), eps);
        prox = frac;
    end

    % Simple rating
    is_good = (snr_db >= snr_good_db) & (stab_std_mm <= stab_good_mm) & (prox <= resprox_good_frac);
    is_ok   = ~is_good & (snr_db >= (snr_good_db-5)) & (stab_std_mm <= (1.5*stab_good_mm));
    rate = 0;
    if is_ok
        rate = 1;
    end
    if is_good
        rate = 2;
    end
    % 0:poor , 1:ok , 2:good

    Q(k,:) = [xs, xe, snr_db, stab_std_mm, prox, rate];
end

Tq = array2table(Q, 'VariableNames', ...
    {'x_start_m','x_end_m','snr_db','bootstrap_std_mm','resonance_prox','rating'});

% Save CSV in THIS script’s folder
writetable(Tq, fullfile(save_folder,'quality.csv'));

quality_output.quality = struct('table', Tq, ...
                    'segment_len_m', segment_length_m, ...
                    'thresholds', struct('snr_good_db',snr_good_db, ...
                                         'stab_good_mm',stab_good_mm, ...
                                         'resprox_good_frac',resprox_good_frac));

% Save the new results to THIS script's output folder
output_path = fullfile(save_folder, 'quality_results.mat');
save(output_path, 'quality_output', '-v7.3');
fprintf('Saved new quality metrics to: %s\n', output_path);

%% --------------------------------- Plots -----------------------------------
if make_plots
    % Fig 1 — Profile with shaded CI
    f1 = mkfig(sprintf('Reconstructed profile with CI [%s]', recon_series), [100 70 1150 520]);
    hold on; grid on; box on; formatAxes();

    fill_ci(x_grid, loZ, hiZ, [0.8 0.88 1.0]);             % shaded CI
    plot(x_grid, zr_base, 'k-', 'LineWidth', 1.1);
    plot(x_grid, mZ,      'b-', 'LineWidth', 1.2);
    if isfield(cfg,'ground') && isfield(cfg.ground,'x') && ~isempty(cfg.ground.x)
        plot(cfg.ground.x, cfg.ground.zr, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
        legend('95% CI','baseline recon','bootstrap mean','ground truth','Location','best');
    else
        legend('95% CI','baseline recon','bootstrap mean','Location','best');
    end
    title('Reconstructed road with 95% confidence band');
    xlabel('Distance (m)'); ylabel('Elevation (m)');
    if export_png
        exportgraphics(f1, fullfile(fig_dir,'38_ci_shaded_profile.png'), 'Resolution', 180);
    end

    % Fig 2 — Segment quality bars
    f2 = mkfig('Segment quality metrics', [130 90 1150 740]);
    tlo2 = tiledlayout(f2, 3, 1, 'TileSpacing','compact','Padding','compact');

    centers = 0.5*(Tq.x_start_m + Tq.x_end_m);
    nexttile(tlo2,1);
    bar(centers, Tq.snr_db, 1.0, 'FaceColor',[0.55 0.7 0.95], 'EdgeColor','none'); grid on; box on; formatAxes();
    yline(snr_good_db, 'k--', 'good threshold');
    ylabel('SNR (dB)'); title('Segment SNR (accel)');

    nexttile(tlo2,2);
    bar(centers, Tq.bootstrap_std_mm, 1.0, 'FaceColor',[0.65 0.85 0.65], 'EdgeColor','none'); grid on; box on; formatAxes();
    yline(stab_good_mm, 'k--', 'good threshold');
    ylabel('Bootstrap std (mm)'); title('Stability across bootstraps');

    nexttile(tlo2,3);
    rcol = [0.80 0.55 0.65];
    plot(centers, Tq.resonance_prox, '-o', 'LineWidth',1.1,'MarkerSize',4,'Color',rcol);
    grid on; box on; formatAxes(); yline(resprox_good_frac, 'k--', 'good threshold');
    xlabel('Distance (m)'); ylabel('Fraction'); title('Resonance proximity (power fraction near wheel-hop)');

    if export_png, exportgraphics(f2, fullfile(fig_dir,'39_segment_quality.png'), 'Resolution', 180); end

    % Fig 3 — Bootstrap std map across x
    f3 = mkfig('Bootstrap std across distance', [160 110 1150 420]);
    plot(x_grid, 1e3*sZ, 'm-', 'LineWidth', 1.1); grid on; box on; formatAxes();
    xlabel('Distance (m)'); ylabel('STD (mm)'); title('Pointwise bootstrap standard deviation');
    if export_png, exportgraphics(f3, fullfile(fig_dir,'40_bootstrap_std_map.png'), 'Resolution', 180); end

    % Fig 4 — Acceleration PSD and noise-estimation band
    f4 = mkfig('Acceleration PSD and noise estimate band', [190 130 980 540]);
    [Pfull, ffull] = pwelch(az_s - mean(az_s,'omitnan'), [], [], [], fs, 'onesided');
    loglog(ffull, Pfull, 'b', 'LineWidth', 1.1); hold on; grid on; box on; formatAxes();
    f_wh_nom = (1/(2*pi))*sqrt(cfg.kt/cfg.mu);
    % Use plot() instead of xline() for broad compatibility
    yl = ylim; plot([f_wh_nom f_wh_nom], yl, '--k', 'LineWidth', 1.0);
    if strcmpi(noise_method,'hf_psd')
        hf_start = hf_band_factor*f_wh_nom;
        plot([hf_start hf_start], yl, ':r', 'LineWidth', 1.0);
        text(hf_start, yl(2), sprintf('  HF start (%.1f×f_{wh})', hf_band_factor), ...
            'VerticalAlignment','top','HorizontalAlignment','left', 'Rotation',90, 'Color',[0.6 0 0]);
    end
    xlabel('Frequency (Hz)'); ylabel('PSD (m^2/s^4/Hz)'); title('Sprung acceleration PSD');
    if export_png, exportgraphics(f4, fullfile(fig_dir,'41_accel_psd_and_noise_est.png'), 'Resolution', 180); end
end

%% =============================== FUNCTIONS ===============================
function zr_x = recon_from_accel_once(t, az, cfg, x_time, x_grid, ...
    hp_fc, hp_or, useSG, win_sec, poly, recon_series, sp)
% End-to-end single-pass reconstruction mirroring the paper method.
% Uses HP drift control on z_s, two LTI steps to z_u then z_r, spatial mapping,
% and (optionally) reuses spatial LP from sp if recon_series='filtered'.

    az = az(:) - mean(az(:),'omitnan');
    t  = t(:); x_time = x_time(:); x_grid = x_grid(:);

    % (1) Integrate acceleration -> velocity -> displacement
    vz_raw = cumtrapz(t, az);
    zs_raw = cumtrapz(t, vz_raw);

    % (2) Drift control high-pass on zs and vz
    if hp_fc > 0
        Wn = clamp01(hp_fc/(cfg.fs/2));
        [b,a] = butter(max(1,hp_or), Wn, 'high');
        zs = filtfilt(b,a, zs_raw);
        vz = filtfilt(b,a, vz_raw);
    else
        zs = zs_raw; vz = vz_raw;
    end

    % (3) Solve C_s z_u̇ + K_s z_u = m_s a_s + c_s z_ṡ + k_s z_s
    cs_eff = max(cfg.cs, 1e-6); % Add protection
    ct_eff = max(cfg.ct, 1e-6); % Add protection

    A1 = -(cfg.ks/cs_eff); B1 = (1/cs_eff); C1 = 1; D1 = 0;
    sys1 = ss(A1,B1,C1,D1);
    zu = lsim(sys1, cfg.ms*az + cs_eff*vz + cfg.ks*zs, t);

    % Smooth derivatives
    vzu = smooth_derivative(zu, t, useSG, win_sec, poly);
    azu = smooth_derivative(vzu, t, useSG, win_sec, poly);

    % (4) Solve tyre/unsprung balance for z_r
    rhs = cfg.mu*azu + (cs_eff+ct_eff)*vzu + (cfg.ks+cfg.kt)*zu - (cs_eff*vz + cfg.ks*zs);
    A2 = -(cfg.kt/ct_eff); B2 = (1/ct_eff); C2 = 1; D2 = 0;
    sys2 = ss(A2,B2,C2,D2);
    zr_t = lsim(sys2, rhs, t);

    % (5) Temporal -> spatial
    zr_x = interp1(x_time, zr_t, x_grid, 'linear', 'extrap'); zr_x(isnan(zr_x)) = 0;

    % (6) Spatial LP if requested (reuse settings from sp)
    if strcmpi(recon_series, 'filtered') && isfield(sp,'nc_spatial') && isfield(sp,'lp_order')
        fs_spatial = 1/max(mean(diff(x_grid),'omitnan'), eps);
        Wn = clamp01(sp.nc_spatial/(fs_spatial/2));
        [b_lp,a_lp] = butter(max(1,sp.lp_order), Wn, 'low');
        zr_x = filtfilt(b_lp, a_lp, zr_x);
    end
end

function v = smooth_derivative(y, t, useSG, win_sec, poly)
% Savitzky–Golay (preferred) or moving-average + gradient fallback.
    y = y(:); t = t(:);
    if numel(y) < 3
        v = zeros(size(y));
        return;
    end
    dt = mean(diff(t),'omitnan');
    if isempty(dt) || dt<=0
        dt = (t(end)-t(1))/max(numel(t)-1,1);
    end

    if useSG && exist('sgolayfilt','file')==2
        w = max(5, 2*floor((win_sec/dt)/2)+1);
        w = min(w, numel(y) - (mod(numel(y),2)==0));
        poly = min(max(poly,2),5);
        y_s = sgolayfilt(y, poly, w);
        v   = gradient(y_s, t);
    else
        w = max(5, 2*floor((win_sec/dt)/2)+1);
        y_s = movmean(y, w, 'Endpoints','shrink');
        v   = gradient(y_s, t);
    end
end

function y = clamp01(y)
    y = max(min(y,0.999), 1e-6);
end

function f = mkfig(name, pos)
    f = figure('Name', name, 'Color', 'w', 'Position', pos);
end

function formatAxes(ax)
    if nargin < 1
        ax = gca;
    end
    ax.LineWidth = 1.0;
    ax.FontName  = 'Calibri';
    ax.FontSize  = 11;
    grid(ax,'on'); box(ax,'on');
end

function fill_ci(x, lo, hi, col)
% Shaded confidence band helper
    x = x(:); lo = lo(:); hi = hi(:);
    [x,idx] = sort(x,'ascend'); lo=lo(idx); hi=hi(idx);
    X = [x; flipud(x)]; Y = [lo; flipud(hi)];
    patch('XData',X,'YData',Y,'FaceColor',col,'EdgeColor','none','FaceAlpha',0.5);
end