
%% 13) Analyze Reconstruction Robustness to Sensor Noise
% -------------------------------------------------------------------------
% PURPOSE:
% This script investigates how robust the road profile reconstruction
% algorithm is to measurement noise from the vehicle's accelerometer. It
% measures the degradation in reconstruction accuracy as the level of
% input noise is increased.
%
% CORE LOGIC (How it works):
% 1.  Loads the configuration file ('recon_cfg.mat') to get the clean input
%     acceleration signal, ground truth profile, and vehicle parameters.
%
% 2.  Defines a series of noise levels to test (e.g., 0.005g, 0.01g, 0.02g,
%     etc.).
%
% 3.  It then loops through each noise level. In each iteration, it
%     generates a random, zero-mean white noise signal with the specified
%     standard deviation and adds it to the clean acceleration signal.
%
% 4.  For each newly created noisy acceleration signal, it re-runs the
%     ENTIRE reconstruction pipeline internally (time-domain inversion and
%     spatial mapping).
%
% 5.  It compares the resulting reconstructed profile to the ground truth
%     and calculates a set of error metrics (e.g., RMSE, correlation).
%
% 6.  After the loop finishes, it saves all the results, showing how the
%     error metrics change as a function of the input noise level.
%
% 7.  It generates plots to visualize this relationship, such as RMSE vs.
%     Noise Level, to clearly show how performance degrades.
%
% INPUTS:
% - '02_Reconstruction Inputs/recon_cfg.mat': The baseline configuration
%   with the clean acceleration signal.
%
% OUTPUTS:
% - '13_Noise Robustness/noise_robustness.mat': Contains the full dataset
%   from the noise sweep.
% - '13_Noise Robustness/noise_robustness.csv': A summary CSV of metrics
%   vs. noise level.
% - Figures visualizing the robustness results.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ========================== HYPERPARAMETERS ==========================
% (1) Input discovery
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};

% (2) Outputs
out_root   = fullfile('00_Outputs', '04_Validation', 'NoiseRobustness');
fig_dir    = fullfile(out_root,'figs');

% (3) Noise sweep in acceleration (σ in g)
noise_levels_g = [0, 0.005, 0.01, 0.02, 0.05, 0.10];  % g-rms
rng_seed       = 42;                                  % reproducibility

% (4) Choose base acceleration to contaminate
%     'clean' -> use clean if available (preferred for isolation)
%     'meas'  -> use cfg.az_s as-is (may already include noise)
base_accel_source = 'clean';   % 'clean' | 'meas'

% (5) Drift control & derivatives (paper defaults)
hp_disp_fc_hz  = 0.05;     % high-pass cutoff for z_s (Hz)
hp_disp_order  = 2;        % IIR order
use_sgolay     = true;     % derivative smoothing method
sgolay_win_sec = 0.15;     % s
sgolay_poly    = 3;

% (6) Spatial low-pass (paper policy)
lp_order          = 4;     % Butterworth order for spatial LP
restore_mean_post = true;  % restore mean after LP

% (7) Tyre damping (small epsilon to avoid singularity if ct=0)
min_ct_Ns_per_m   = 1.0;   % Ns/m

% (8) Plots & export
overlay_levels_g = [0, 0.02, 0.10];  % which noise levels to overlay vs truth
make_plots       = true;
export_png       = true;

%% ------------------------- Load configuration -------------------------
cfg_path = first_existing(cfg_candidates);
assert(~isempty(cfg_path), 'recon_cfg.mat not found in expected locations.');
C = load(cfg_path); cfg = C.cfg;

if ~exist(out_root,'dir'), mkdir(out_root); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Ground truth required
assert(isfield(cfg,'ground') && isfield(cfg.ground,'x') && ~isempty(cfg.ground.x), ...
    'Ground-truth profile not available in cfg.ground. Run your simulator to include it.');

% Time/grid/vehicle
t   = cfg.t(:);             fs = cfg.fs;   dt = cfg.dt;
V   = cfg.V;                x_t = cfg.x_t(:);
x_true  = cfg.ground.x(:);  zr_true = cfg.ground.zr(:);

% Base acceleration to contaminate
% NOTE: For this to be the "clean" signal, ensure 'accel_channel' was set to 'clean' in recon_init_and_load.m
az_base = cfg.az_s(:);

% Parameters (Golden Car / user cfg)
P.ms = cfg.ms; P.mu = cfg.mu; P.ks = cfg.ks; P.cs = cfg.cs; P.kt = cfg.kt;
P.ct = max(cfg.ct, min_ct_Ns_per_m);   % ensure > 0 for numerical stability

% Spatial evaluation grid: use ground-truth x for metrics
x_grid = x_true;
dx_grid = mean(diff(x_grid),'omitnan');
fs_spatial = 1 / max(dx_grid, eps);

%% ------------------------- Sweep over noise levels -------------------------
rng(rng_seed);
g0 = 9.80665;

N = numel(noise_levels_g);
metrics = struct('sigma_g',[],'sigma_mps2',[],'snr_db',[], ...
                 'rmse',[],'mae',[],'bias',[],'r',[],'rse_all',[]);
profiles = cell(N,1);   % store reconstructed profiles per level

% Keep PSDs of accel for the lowest and highest noise, for plotting
accel_psd = struct('f',[],'Pclean',[],'Pnoisy_lo',[],'Pnoisy_hi',[]);

% Clean acceleration for SNR reference
az_clean_ref = az_base(:);  % treat base as "clean" reference for SNR (if base already noisy, SNR is relative)

% PSD of the clean reference (for Fig 37)
[accel_psd.Pclean, accel_psd.f] = pwelch_centered(az_clean_ref, fs);

% Reconstruction & metrics loop
for i = 1:N
    sigma_g = noise_levels_g(i);
    sigma_mps2 = sigma_g * g0;

    % Generate and add noise (white Gaussian)
    noise = sigma_mps2 * randn(size(az_base));
    az_noisy = az_base + noise;

    % SNR estimate (time-domain rms)
    snr_db = 10*log10( rms_ex(az_clean_ref)^2 / max(rms_ex(noise)^2, eps) );

    % Reconstruct once
    rec = reconstruct_once(t, az_noisy, fs, V, x_t, x_grid, P, ...
                           hp_disp_fc_hz, hp_disp_order, use_sgolay, sgolay_win_sec, sgolay_poly, ...
                           lp_order, restore_mean_post);

    profiles{i} = struct('sigma_g',sigma_g,'x',x_grid,'zr',rec.zr_x,'zr_raw',rec.zr_x_raw,'nc',rec.nc);

    % Metrics vs truth
    M = compare_metrics(x_grid, rec.zr_x, x_true, zr_true, fs_spatial, true);

    % Record
    metrics(i).sigma_g     = sigma_g;
    metrics(i).sigma_mps2  = sigma_mps2;
    metrics(i).snr_db      = snr_db;
    metrics(i).rmse        = M.pointwise.rmse_m;
    metrics(i).mae         = M.pointwise.mae_m;
    metrics(i).bias        = M.pointwise.bias_m;
    metrics(i).r           = M.pointwise.pearson_r;
    metrics(i).rse_all     = M.spectral.rse_all;

    % Store accel PSDs for extremes
    if i==1
        [accel_psd.Pnoisy_lo, ~] = pwelch_centered(az_noisy, fs);
    end
    if i==N
        [accel_psd.Pnoisy_hi, ~] = pwelch_centered(az_noisy, fs);
    end
end

%% ------------------------- Save results -------------------------
R = struct();
R.settings.noise_levels_g   = noise_levels_g;
R.settings.rng_seed         = rng_seed;
R.settings.base_accel_src   = base_accel_source;
R.settings.hp_fc_hz         = hp_disp_fc_hz;
R.settings.hp_order         = hp_disp_order;
R.settings.use_sgolay       = use_sgolay;
R.settings.sgolay_win_sec   = sgolay_win_sec;
R.settings.sgolay_poly      = sgolay_poly;
R.settings.lp_order         = lp_order;
R.settings.restore_mean     = restore_mean_post;
R.settings.min_ct_Ns_per_m  = min_ct_Ns_per_m;

R.params   = P;
R.profiles = profiles;
R.metrics  = metrics;
R.accel_psd = accel_psd;

save(fullfile(out_root,'noise_robustness.mat'), 'R','-v7.3');
fprintf('Saved: %s\n', fullfile(out_root,'noise_robustness.mat'));

% CSV summary
csv_path = fullfile(out_root,'noise_robustness.csv');
fid = fopen(csv_path,'w');
assert(fid>0, 'Cannot open %s for writing', csv_path);
fprintf(fid,'sigma_g,sigma_mps2,snr_db,rmse_m,mae_m,bias_m,r,rse_all\n');
for i=1:N
    fprintf(fid,'%.6f,%.9e,%.3f,%.9e,%.9e,%.9e,%.6f,%.6f\n', ...
        metrics(i).sigma_g, metrics(i).sigma_mps2, metrics(i).snr_db, ...
        metrics(i).rmse, metrics(i).mae, metrics(i).bias, metrics(i).r, metrics(i).rse_all);
end
fclose(fid);
fprintf('Saved: %s\n', csv_path);

%% ------------------------- Plots -------------------------
if make_plots
    % --- Fig 1: Metrics vs noise level ---
    f1 = mkfig('Noise robustness — metrics vs noise', [100 70 1160 760]);
    tlo1 = tiledlayout(f1,3,1,'TileSpacing','compact','Padding','compact');

    sig = [metrics.sigma_g];  % g
    nexttile(tlo1,1);
    plot(sig, [metrics.rmse], '-o','LineWidth',1.3,'MarkerSize',4); grid on; box on; formatAxes();
    title('RMSE vs noise level'); xlabel('Noise \sigma_a (g)'); ylabel('RMSE of profile (m)');

    nexttile(tlo1,2);
    plot(sig, [metrics.r], '-o','LineWidth',1.3,'MarkerSize',4,'Color',[0.25 0.45 0.85]); grid on; box on; formatAxes();
    title('Pearson r vs noise level'); xlabel('Noise \sigma_a (g)'); ylabel('Correlation r');

    nexttile(tlo1,3);
    plot(sig, [metrics.rse_all], '-o','LineWidth',1.3,'MarkerSize',4,'Color',[0.60 0.20 0.70]); grid on; box on; formatAxes();
    title('Spectral RSE (overall) vs noise level'); xlabel('Noise \sigma_a (g)'); ylabel('RSE');

    if export_png, exportgraphics(f1, fullfile(fig_dir,'35_noise_metrics_vs_level.png'), 'Resolution',180); end

    % --- Fig 2: Truth vs recon overlays for selected noise levels ---
    f2 = mkfig('Noise robustness — overlays', [140 100 1180 820]);
    tlo2 = tiledlayout(f2, numel(overlay_levels_g), 1, 'TileSpacing','compact','Padding','compact');

    for k=1:numel(overlay_levels_g)
        % Find nearest in noise_levels_g
        [~, idx] = min(abs([metrics.sigma_g] - overlay_levels_g(k)));
        Pk = profiles{idx};
        nexttile(tlo2,k);
        plot(x_true, zr_true, 'k', 'LineWidth', 1.2); hold on;
        plot(Pk.x, Pk.zr, 'b', 'LineWidth', 1.1);
        grid on; box on; formatAxes();
        title(sprintf('Overlay at \\sigma_a = %.3f g  (RMSE = %.2f mm, r = %.3f, RSE = %.3f)', ...
              metrics(idx).sigma_g, 1e3*metrics(idx).rmse, metrics(idx).r, metrics(idx).rse_all));
        xlabel('Distance (m)'); ylabel('Elevation (m)');
        legend('Truth','Recon','Location','best');
    end
    if export_png, exportgraphics(f2, fullfile(fig_dir,'36_overlay_selected_levels.png'), 'Resolution',180); end

    % --- Fig 3: Input acceleration PSD (clean vs extremes) ---
    f3 = mkfig('Noise robustness — input acceleration PSD', [170 130 1000 640]);
    loglog(accel_psd.f, accel_psd.Pclean, 'k','LineWidth',1.2); hold on;
    if ~isempty(accel_psd.Pnoisy_lo)
        loglog(accel_psd.f, accel_psd.Pnoisy_lo, 'b','LineWidth',1.1);
    end
    if ~isempty(accel_psd.Pnoisy_hi)
        loglog(accel_psd.f, accel_psd.Pnoisy_hi, 'm','LineWidth',1.1);
    end
    grid on; box on; formatAxes();
    title('PSD of sprung acceleration a_z^s (reference vs noisy)'); xlabel('Frequency (Hz)'); ylabel('PSD (m^2/s^4/Hz)');
    legend('reference','low-noise case','high-noise case','Location','southwest');
    if export_png, exportgraphics(f3, fullfile(fig_dir,'37_input_accel_psd.png'), 'Resolution',180); end
end

%% =============================== FUNCTIONS ===============================
function p = first_existing(candidates)
    p = '';
    for i = 1:numel(candidates)
        if exist(candidates{i},'file')==2, p = candidates{i}; return; end
    end
end

function rec = reconstruct_once(t, azs, fs, V, x_t, x_grid, P, ...
                                hp_fc, hp_or, useSG, win_sec, poly, ...
                                lp_order, restore_mean)
    % End-to-end reconstruction for a given acceleration record and parameters.
    azs = azs(:) - mean(azs,'omitnan');
    t   = t(:); x_t = x_t(:); x_grid = x_grid(:);

    % Integrate to z_s, drift-control high-pass (paper: fc≈0.05 Hz)
    vzs_raw = cumtrapz(t, azs);
    zss_raw = cumtrapz(t, vzs_raw);
    [zss, vzs] = hp_disp(zss_raw, vzs_raw, fs, hp_fc, hp_or);

    % Internal sprung force f_s = m_s a_s + c_s z_ṡ + k_s z_s
    fs_in = P.ms*azs + P.cs*vzs + P.ks*zss;

    % Step 1: C_s z_u̇ + K_s z_u = f_s  → z_u via first-order LTI
    cs_eff = max(P.cs, 1e-6); % Add protection
    ct_eff = max(P.ct, 1e-6); % Add protection

    A1 = -(P.ks/cs_eff); B1 = (1/cs_eff); C1 = 1; D1 = 0;
    sys1 = ss(A1,B1,C1,D1);
    zu   = lsim(sys1, fs_in, t);

    % Derivatives of z_u with smoothing
    vzu = smooth_derivative(zu, t, useSG, win_sec, poly);
    azu = smooth_derivative(vzu, t, useSG, win_sec, poly);

    % Step 2: C_t z_ṙ + K_t z_r = μ a_u + (C_s+C_t) z_u̇ + (K_s+K_t) z_u - (C_s z_ṡ + K_s z_s)
    rhs = P.mu*azu + (cs_eff+ct_eff)*vzu + (P.ks+P.kt)*zu - (P.cs*vzs + P.ks*zss);

    A2 = -(P.kt/ct_eff); B2 = (1/ct_eff); C2 = 1; D2 = 0;
    sys2 = ss(A2,B2,C2,D2);
    zr_t = lsim(sys2, rhs, t);

    % Temporal → spatial and low-pass with n_c = f_wh/(3 V)
    zr_x_raw = interp1(x_t, zr_t, x_grid, 'linear', 'extrap'); zr_x_raw(isnan(zr_x_raw)) = 0;

    dx         = mean(diff(x_grid),'omitnan');
    fs_spatial = 1/max(dx,eps);
    nyq_sp     = fs_spatial/2;

    f_wh = (1/(2*pi))*sqrt(P.kt/P.mu);     % Hz
    nc   = f_wh / 3 / max(V,eps);          % cycles/m
    Wn   = clamp01(nc/nyq_sp);
    [b_lp,a_lp] = butter(max(1,lp_order), Wn, 'low');
    zr_x = filtfilt(b_lp,a_lp, zr_x_raw);

    if restore_mean
        zr_x = zr_x + (mean(zr_x_raw,'omitnan') - mean(zr_x,'omitnan'));
    end

    rec = struct('zr_x',zr_x(:), 'zr_x_raw',zr_x_raw(:), 'nc',nc, 'fs_spatial',fs_spatial);
end

function [M] = compare_metrics(x_rec, zr_rec, x_true, zr_true, fs_spatial, remove_dc_for_psd)
    % Interpolate reconstruction onto truth grid and compute metrics & spectral RSE.
    xmin = max(min(x_true), min(x_rec));
    xmax = min(max(x_true), max(x_rec));
    mT = (x_true>=xmin) & (x_true<=xmax);
    mR = (x_rec >=xmin) & (x_rec <=xmax);
    xT = x_true(mT); yT = zr_true(mT);
    xR = x_rec(mR);  yR = zr_rec(mR);

    y_rec_on_true = safe_interp(xR, yR, xT);
    e = y_rec_on_true - yT;

    M.pointwise.rmse_m = sqrt(mean(e.^2, 'omitnan'));
    M.pointwise.mae_m  = mean(abs(e), 'omitnan');
    M.pointwise.bias_m = mean(e, 'omitnan');
    M.pointwise.pearson_r = robust_corr(yT, y_rec_on_true);

    % Welch PSD (spatial), overall RSE
    [Ptrue, Ntrue] = welch_psd_spatial(yT, fs_spatial, [], 0.5, 0.5, remove_dc_for_psd);
    [Prec,  Nrec ] = welch_psd_spatial(y_rec_on_true, fs_spatial, [], 0.5, 0.5, remove_dc_for_psd);
    n_common = linspace(0, min([Ntrue(end), Nrec(end), fs_spatial/2]), max(numel(Ntrue),numel(Nrec)))';
    Ptrue_c = safe_interp(Ntrue, Ptrue, n_common);
    Prec_c  = safe_interp(Nrec,  Prec,  n_common);
    M.spectral.rse_all = sqrt( sum((Prec_c - Ptrue_c).^2) / (sum(Ptrue_c.^2) + eps) );
end

function [Pxx, N, cfg_out] = welch_psd_spatial(x, fs_spatial, nfft_user, win_frac, ovlp_frac, remove_dc)
    x = x(:);
    if remove_dc, x = x - mean(x,'omitnan'); end
    Npts = numel(x);
    if isempty(nfft_user)
        nfft = max(1024, 2^nextpow2(min(8192, Npts)));
    else
        nfft = max(256, nfft_user);
    end
    wlen  = max(64, round(win_frac * nfft));   wlen  = min(wlen, Npts);
    if mod(wlen,2)==0, wlen = wlen-1; end
    nover = max(0, round(ovlp_frac * wlen));
    win   = hamming(wlen);
    [Pxx, N] = pwelch(x, win, nover, nfft, fs_spatial, 'onesided');
    cfg_out = struct('nfft',nfft,'wlen',wlen,'nover',nover,'remove_dc',logical(remove_dc));
end

function [P,f] = pwelch_centered(x, fs)
    x = x(:) - mean(x,'omitnan');
    N = max(1024, 2^nextpow2(min(16384, numel(x))));
    w = hamming(round(0.5*N));
    nover = round(0.5*length(w));
    [P,f] = pwelch(x, w, nover, N, fs, 'onesided');
end

function [z_hp, vz_hp] = hp_disp(z_raw, vz_raw, fs, fc, order)
    if isempty(fc) || fc<=0
        z_hp = z_raw; vz_hp = vz_raw; return;
    end
    Wn = clamp01(fc/(fs/2));
    [b,a] = butter(max(1,order), Wn, 'high');
    z_hp  = filtfilt(b,a, z_raw);
    vz_hp = filtfilt(b,a, vz_raw);
end

function v = smooth_derivative(y, t, useSG, win_sec, poly)
    y = y(:); t = t(:);
    if numel(y) < 3, v = zeros(size(y)); return; end
    dt = mean(diff(t),'omitnan');
    if isempty(dt) || dt<=0, dt = (t(end)-t(1))/max(numel(t)-1,1); end
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

function yi = safe_interp(x, y, xi)
    x = x(:); y = y(:); xi = xi(:);
    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    if numel(x) < 2
        yi = nan(size(xi));
        return;
    end
    [x, idx] = sort(x); y = y(idx);
    yi = interp1(x, y, xi, 'linear', 'extrap');
end

function r = robust_corr(a, b)
    a = a(:); b = b(:);
    m = isfinite(a) & isfinite(b);
    a = a(m); b = b(m);
    if numel(a) < 3
        r = NaN;
        return;
    end
    C = corrcoef(a, b);
    if all(size(C) >= 2) && isfinite(C(1,2))
        r = C(1,2);
    else
        r = NaN;
    end
end

function y = clamp01(y)
    y = max(min(y,0.999), 1e-6);
end

function r = rms_ex(x)
    x = x(:);
    r = sqrt(mean(x.^2, 'omitnan'));
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