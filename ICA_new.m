
%% 3C) Reconstruct Road Profile Using Independent Component Analysis (ICA)
% -------------------------------------------------------------------------
% PURPOSE:
% Reconstructs the road profile z_r(t) (and spatial z_r(x)) from ONLY the
% sprung vertical acceleration a_s(t) using an ICA-based separation workflow,
% following the logic described in:
%   "Road Profile Identification Using Estimation Techniques: Comparison
%    Between Independent Component Analysis and Kalman Filter" (JTAM, 2019).
%
% CORE LOGIC (paper-aligned):
% 1. Load: Reads `recon_cfg.mat` (sprung accel a_s, speed V, time base, params).
% 2. Form Observation: O(t) := a_s(t) (paper uses sprung acceleration as observation).
% 3. Preprocess: Center + Whiten the observation (required by ICA).
% 4. ICA Separation: Estimate separating matrix/vector W by maximizing
%    non-Gaussianity (kurtosis) using a deflation-style extraction.
% 5. Source -> Road: Treat the extracted independent component as the road
%    profile estimate (up to an unknown scale/sign factor; ICA ambiguity).
% 6. Postprocess: Optional detrend/band-limit, then map z_r(t) -> z_r(x).
% 7. Output: Saves standard spatial struct `sp` to
%    '03C_ICA/recon_spatial.mat' and optionally publishes to shared folder.
%
% NOTES / PRACTICALITIES:
% - With a single measured channel, ICA requires additional structure to be
%   non-trivial; this implementation uses a time-delay embedding of a_s(t)
%   to create a multichannel observation matrix before ICA, while preserving
%   the paper’s core steps (centering, whitening, kurtosis-based W, deflation).
% - The estimated road profile is recovered up to a scale factor. The code
%   provides optional amplitude calibration methods (disabled by default).
% -------------------------------------------------------------------------

close all; clc;
fprintf('\n=== PATH C: ICA-based reconstruction ===\n');

%% ---------------------------- HYPERPARAMETERS -----------------------------
% --- Input / observation formation ---
obs_signal_name       = 'sprung_accel'; % for figure titles/metadata only

% --- Time-delay embedding (to make ICA feasible with 1 sensor) ---
use_delay_embedding   = true;   % recommended: true for single-channel input
embed_dim             = 12;     % embedding dimension M (channels) [typ. 8–30]
embed_delay_samples   = 2;      % delay tau in samples between embedded channels
embed_crop_policy     = 'valid';% {'valid','same'} valid -> shorter, consistent

% --- ICA optimization (kurtosis-based, deflation) ---
ica_max_iter          = 1000;   % max iterations for each component
ica_tol               = 1e-6;   % convergence threshold on w change
ica_n_components      = 1;      % recover 1 component interpreted as road input
ica_contrast          = 'kurtosis'; % fixed to paper logic
ica_nonlin            = 'pow3'; % y^3 nonlinearity for kurtosis fixed-point

% --- ICA solver selection ---
ica_solver            = 'fastica_fp';  % {'robustica_os','fastica_fp'}
robustica_target_sign = 0;              % 0 = maximize |kurtosis| (paper), +1 or -1 to target sign
robustica_normalize_g = true;           % recommended by RobustICA for conditioning

ica_orthogonalize     = true;   % deflation orthogonalization between components
ica_seed              = 1;      % RNG seed for reproducibility
ica_init              = 'random_unit'; % {'random_unit','first_pc'}

% --- Post-processing on recovered source (optional but useful) ---
do_detrend_source     = true;   % remove mean/linear trend in estimated road
do_bandlimit_source   = true;   % band-limit in frequency domain (recommended)
f_hp                  = 0.05;   % high-pass cutoff [Hz] (remove drift)
f_lp                  = 80;     % low-pass cutoff [Hz] (remove HF noise)
bandlimit_order       = 4;      % Butterworth order (time-domain filter option)
bandlimit_method      = 'fft';  % {'fft','butter'} (fft preserves phase nicely)

% --- Amplitude / scale calibration (ICA scale ambiguity) ---
calibrate_scale       = 'match_rms'; % {'none','match_rms','match_psd','match_truth'}
target_rms_m          = 0.0028;  % [m] only used if calibrate_scale='match_rms'
% NOTE: 'match_truth' requires cfg.ground.zr to exist.

% --- Spatial mapping ---
spatial_grid_policy   = 'ground_truth'; % {'ground_truth','uniform'}

% --- Fig C4: binned RMSE vs distance ---
rmse_bin_m            = 1;     % bin size in meters (e.g., 10 m)
rmse_bin_min_points   = 5;      % minimum samples required in a bin to compute RMSE
rmse_bin_plot_style   = 'line'; % {'stairs','bar','line'}
rmse_bin_use_robust_x = true;   % true -> use x_common for binning (aligned grid)

% --- Output / publishing ---
publish_to_shared     = true;
out_root_dir          = '00_Outputs';
out_dir_pathC         = fullfile(out_root_dir, '03_ReconstructionCore', 'PathC_ICA');

% --- Plot controls ---
make_plots            = true;
save_plots            = true;     % write PNG + FIG under out_dir_pathC/figs
fig_format            = 'png';    % {'png','tiff','pdf'}
fig_dpi               = 300;
use_latex_labels      = false;

% Plot style
baseFont              = 11;
lineLW                = 1.6;
useMinorGrid          = true;

% --- Diagnostics ---
print_ica_diagnostics = true;    % prints iterations, kurtosis, convergence

%% ----------------------------- OUTPUT FOLDERS -----------------------------
fig_dir = fullfile(out_dir_pathC,'figs');
if ~exist(out_dir_pathC,'dir')
    mkdir(out_dir_pathC);
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
assert(~isempty(cfg_path), 'PATH C: recon_cfg.mat not found. Run prepare_reconstruction_inputs.m first.');

S = load(cfg_path);
assert(isfield(S,'cfg'), 'PATH C: recon_cfg.mat missing cfg struct.');
cfg = S.cfg;

% --- Time base and sampling ---
t  = cfg.t(:);
fs = cfg.fs;

% --- Vehicle speed ---
assert(isfield(cfg,'V') && ~isempty(cfg.V), 'PATH C: cfg.V (speed) missing.');
V = cfg.V;

% --- Locate sprung acceleration (observation O(t) := a_s(t)) ---
accel_field_order = {'az_s','a_z_s','a_s','azs','a_sprung','sprung_accel','acc_sprung'};
az_s = [];
for f = accel_field_order
    if isfield(cfg, f{1}) && ~isempty(cfg.(f{1}))
        az_s = cfg.(f{1})(:);
        break;
    end
end
assert(~isempty(az_s), ...
    'PATH C: Could not locate sprung acceleration in cfg (tried %s).', ...
    strjoin(accel_field_order,', '));

% --- Quarter-car parameters (kept for metadata consistency; ICA itself does not need them) ---
reqParams = {'ms','mu','ks','cs','kt'};
for p = reqParams
    assert(isfield(cfg,p{1}) && ~isempty(cfg.(p{1})), 'PATH C: cfg.%s missing.', p{1});
end
ms = cfg.ms;
mu = cfg.mu;
ks = cfg.ks;
cs = cfg.cs;
kt = cfg.kt;

% Optional tire damping if present
ct = 0;
if isfield(cfg,'ct') && ~isempty(cfg.ct)
    ct = cfg.ct;
end

% --- Distance axis for mapping ---
if isfield(cfg,'x_t') && ~isempty(cfg.x_t)
    x_t = cfg.x_t(:);
else
    x_t = (t - t(1)) * V;
end

% --- Basic checks ---
N = numel(t);
assert(numel(az_s) == N, 'PATH C: az_s length (%d) != t length (%d).', numel(az_s), N);

fprintf('Inputs: N=%d, fs=%.3f Hz, V=%.3f m/s\n', N, fs, V);

% Observation (paper logic)
O = az_s;   % O(t) := a_s(t)

%% -------------------- 1) Build ICA Observation Matrix ---------------------
% Paper observation: O(t) := a_s(t). For single-channel ICA to be non-trivial,
% we create a multichannel observation via time-delay embedding:
%   Xobs = [a(t), a(t-tau), a(t-2tau), ...]^T
%
% Dimensions:
%   Xobs : [M x N_eff]   (M = embed_dim)

rng(ica_seed);

if use_delay_embedding
    [Xobs, t_eff, x_eff] = build_delay_embedding(O, t, x_t, embed_dim, embed_delay_samples, embed_crop_policy);
else
    % Degenerate 1xN case (ICA cannot separate beyond scale/sign)
    Xobs  = O(:).';      % [1 x N]
    t_eff = t(:);
    x_eff = x_t(:);
end

[M, N_eff] = size(Xobs);
fprintf('ICA observation matrix: M=%d channels, N_eff=%d samples\n', M, N_eff);

assert(N_eff >= max(100, 5*M), ...
    'PATH C: Not enough samples after embedding (N_eff=%d). Reduce embed_dim or embed_delay_samples.', N_eff);

%% -------------------- 2) Center + Whiten (ICA prerequisites) --------------
% Center each channel: remove mean
Xc = Xobs - mean(Xobs, 2);

% Whitening: Z = V^{-1/2} E^T Xc  such that cov(Z) = I
% Returns:
%   Z        : whitened data [M x N_eff]
%   Wwhite   : whitening matrix  [M x M]
%   Dewhite  : dewhitening matrix [M x M]  (inverse map)
[Z, Wwhite, Dewhite, covX] = whiten_channels(Xc);

% Sanity check: covariance should be ~I
if print_ica_diagnostics
    covZ = (Z*Z.') / max(N_eff-1,1);
    offdiag_rms = sqrt(mean((covZ(~eye(size(covZ))).^2), 'omitnan'));
    diag_dev_rms = sqrt(mean((diag(covZ)-1).^2, 'omitnan'));
    fprintf('Whitening check: offdiag RMS=%.3e, diag dev RMS=%.3e\n', offdiag_rms, diag_dev_rms);
end

% Optional: clip extreme outliers (numerical stability, off by default)
% z_clip = 0; % e.g., 6 -> clip to +/-6 std
% if z_clip > 0
%     Z = max(min(Z, z_clip), -z_clip);
% end

%% -------------------- 3) Initialize ICA / Deflation Setup -----------------
% We will estimate W_ica (n_components x M) such that:
%   Y = W_ica * Z
% Each row of W_ica yields one independent component y_k (unit variance).
nComp = max(1, min(ica_n_components, M));
W_ica = zeros(nComp, M);
ica_info = struct();
ica_info.iter = zeros(nComp,1);
ica_info.kurt = zeros(nComp,1);
ica_info.converged = false(nComp,1);

% Initial direction(s)
switch lower(ica_init)
    case 'random_unit'
        w0 = randn(M,1); w0 = w0 / norm(w0);
    case 'first_pc'
        % Use first eigenvector of cov(Z) ~ I -> fall back to random
        w0 = randn(M,1); w0 = w0 / norm(w0);
    otherwise
        error('PATH C: Unknown ica_init: %s', ica_init);
end

% Placeholder for extracted components (time series)
Y = zeros(nComp, N_eff);

for kComp = 1:nComp

    % --- initialize w for this component ---
    if kComp == 1
        w = w0;
    else
        w = randn(M,1); w = w / norm(w);
        if ica_orthogonalize && kComp > 1
            w = deflate_orth(w, W_ica(1:kComp-1,:));
            w = w / max(norm(w), eps);
        end
    end

    converged = false;
    it = 0;

    for it = 1:ica_max_iter
        w_old = w;

        % Project current output
        y = w.' * Z;   % [1 x N_eff]

        switch lower(ica_solver)
            case 'robustica_os'
                % ---------- RobustICA (Zarzoso & Comon): exact line search on kurtosis ----------
                % Gradient direction (real-valued kurtosis case; matches your current data)
                g = (Z * (y.^3).') / N_eff - 3*w;  % proportional to ∇ kurtosis under whitening

                if robustica_normalize_g
                    g = g / max(norm(g), eps);
                end

                % Deflation orthogonalization on direction (keeps search in orthogonal subspace)
                if ica_orthogonalize && kComp > 1
                    g = deflate_orth(g, W_ica(1:kComp-1,:));
                    g = g / max(norm(g), eps);
                end

                % Compute optimal step size mu by quartic polynomial roots (RobustICA)
                mu_opt = robustica_opt_mu_kurtosis_real(w, g, Z, robustica_target_sign);

                % Update
                w = w + mu_opt * g;

            case 'fastica_fp'
                % ---------- Your existing FastICA-like fixed-point ----------
                w = (Z * (y.^3).') / N_eff - 3*w;

            otherwise
                error('Unknown ica_solver: %s', ica_solver);
        end

        % Deflation orthogonalization (paper mentions deflation)
        if ica_orthogonalize && kComp > 1
            w = deflate_orth(w, W_ica(1:kComp-1,:));
        end

        % Normalize
        w = w / max(norm(w), eps);

        % Convergence: direction stabilized (sign-invariant)
        if 1 - abs(w.' * w_old) < ica_tol
            converged = true;
            break;
        end
    end

    % Store separating vector
    W_ica(kComp,:) = w(:).';

    % Extract component
    y = w.' * Z;                          % [1 x N_eff]
    Y(kComp,:) = y;

    % Kurtosis diagnostic (excess kurtosis)
    y0 = y - mean(y);
    m2 = mean(y0.^2);
    m4 = mean(y0.^4);
    if m2 > 0
        kurt_excess = m4/(m2^2) - 3;
    else
        kurt_excess = NaN;
    end

    ica_info.iter(kComp) = it;
    ica_info.kurt(kComp) = kurt_excess;
    ica_info.converged(kComp) = converged;

    if print_ica_diagnostics
        fprintf('ICA comp %d/%d: iter=%d, converged=%d, excess kurtosis=%.4g\n', ...
            kComp, nComp, it, converged, kurt_excess);
    end
end

W_sep = W_ica * Wwhite;    % separation matrix acting on centered observation Xc
Y     = W_sep * Xc;        % Y = W * O_centered  (paper-style)
zr_t_ica = Y(1,:).';

%% -------------------- 5) Source -> Road + Postprocessing ------------------
% Paper interpretation: recovered independent component corresponds to road
% profile (up to scale/sign). We use the first component as z_r(t).

zr_t_ica = Y(1,:).';      % column, length N_eff
t_ica    = t_eff(:);
x_ica    = x_eff(:);

% Optional detrend (remove mean/linear drift)
if do_detrend_source
    zr_t_ica = detrend(zr_t_ica, 1);  % linear detrend
end

% Optional band-limiting (recommended for stability/interpretability)
if do_bandlimit_source
    switch lower(bandlimit_method)
        case 'fft'
            zr_t_ica = fft_bandlimit(zr_t_ica, fs, f_hp, f_lp);
        case 'butter'
            zr_t_ica = butter_bandlimit(zr_t_ica, fs, f_hp, f_lp, bandlimit_order);
        otherwise
            error('PATH C: Unknown bandlimit_method: %s', bandlimit_method);
    end
end

% Scale calibration (ICA ambiguity)
scale_alpha = 1.0;
switch lower(calibrate_scale)
    case 'none'
        % do nothing
    case 'match_rms'
        r_rms = sqrt(mean(zr_t_ica.^2,'omitnan'));
        if isfinite(r_rms) && r_rms > 0
            scale_alpha = target_rms_m / r_rms;
        end
    case 'match_truth'
        if isfield(cfg,'ground') && isfield(cfg.ground,'t') && isfield(cfg.ground,'zr')
            tg = cfg.ground.t(:);
            zg = cfg.ground.zr(:);
            zg_i = interp1(tg, zg, t_ica, 'linear', 'extrap');
            r_rms = sqrt(mean(zr_t_ica.^2,'omitnan'));
            g_rms = sqrt(mean(zg_i.^2,'omitnan'));
            if isfinite(r_rms) && r_rms > 0 && isfinite(g_rms)
                scale_alpha = g_rms / r_rms;
                % Optionally also fix sign to maximize correlation
                if std(zr_t_ica) > 0 && std(zg_i) > 0
                    if corr(zr_t_ica, zg_i, 'Rows','complete') < 0
                        scale_alpha = -scale_alpha;
                    end
                end
            end
        else
            warning('PATH C: match_truth requested but cfg.ground.t/zr missing. Using scale_alpha=1.');
        end
    case 'match_psd'
        warning('PATH C: match_psd not implemented in this section. Using scale_alpha=1.');
    otherwise
        error('PATH C: Unknown calibrate_scale: %s', calibrate_scale);
end

zr_t_ica = scale_alpha * zr_t_ica;

%% --------------------- 6) Map to Spatial + Save Output --------------------
% Map reconstructed road z_r(t) onto distance x, then resample to chosen grid.
% Output artifact matches your standard `sp` struct and saves to Path C folder.

% --- Choose spatial grid ---
switch lower(spatial_grid_policy)
    case 'ground_truth'
        if isfield(cfg,'ground') && isfield(cfg.ground,'x') && ~isempty(cfg.ground.x)
            x_grid = cfg.ground.x(:);
        else
            warning('PATH C: ground truth grid missing; falling back to UNIFORM.');
            x_grid = (t_ica - t_ica(1)) * V;
        end
    case 'uniform'
        x_grid = (t_ica - t_ica(1)) * V;
    otherwise
        error('PATH C: Unknown spatial_grid_policy: %s', spatial_grid_policy);
end

% --- Interpolate from time-based x_ica to x_grid ---
zr_x_raw  = interp1(x_ica, zr_t_ica, x_grid, 'linear', 'extrap');
zr_x_raw(~isfinite(zr_x_raw)) = 0;

% (Keep a "filtered" field for consistency; ICA output is already postprocessed)
zr_x_filt = zr_x_raw;

dx         = mean(diff(x_grid),'omitnan');
fs_spatial = 1 / max(dx, eps);

% --- Build spatial struct (same schema as Path B) ---
sp = struct();
sp.x              = x_grid;
sp.zr_x_raw       = zr_x_raw;
sp.zr_x_filt      = zr_x_filt;
sp.dx             = dx;
sp.fs_spatial     = fs_spatial;
sp.nc_spatial     = NaN;
sp.lp_order       = NaN;
sp.filter_applied = do_bandlimit_source;
sp.meta = struct( ...
    'recon_method',           'Path C: ICA (kurtosis, deflation)', ...
    'fs',                     fs, ...
    'V',                      V, ...
    'spatial_grid_policy',    spatial_grid_policy, ...
    'use_delay_embedding',    use_delay_embedding, ...
    'embed_dim',              embed_dim, ...
    'embed_delay_samples',    embed_delay_samples, ...
    'embed_crop_policy',      embed_crop_policy, ...
    'ica_max_iter',           ica_max_iter, ...
    'ica_tol',                ica_tol, ...
    'ica_n_components',       nComp, ...
    'ica_contrast',           ica_contrast, ...
    'ica_nonlin',             ica_nonlin, ...
    'ica_orthogonalize',      ica_orthogonalize, ...
    'ica_seed',               ica_seed, ...
    'do_detrend_source',      do_detrend_source, ...
    'do_bandlimit_source',    do_bandlimit_source, ...
    'bandlimit_method',       bandlimit_method, ...
    'f_hp',                   f_hp, ...
    'f_lp',                   f_lp, ...
    'scale_alpha',            scale_alpha);

% Save Path C artifact
out_mat_pathC = fullfile(out_dir_pathC,'recon_spatial.mat');
save(out_mat_pathC, 'sp','-v7.3');
fprintf('PATH C: saved %s\n', out_mat_pathC);

% Optional publish for downstream compatibility
if publish_to_shared
    share_dir = fullfile(out_root_dir, '03_ReconstructionCore', 'SpatialMapping');
    if ~exist(share_dir,'dir')
        mkdir(share_dir);
    end
    copyfile(out_mat_pathC, fullfile(share_dir,'recon_spatial.mat'));
    fprintf('PATH C: also published to %s\n', share_dir);
end

%% ------------------------------ 7) PLOTS ----------------------------------
if make_plots
    % ---------- Global style ----------
    if use_latex_labels
        set(0,'defaultTextInterpreter','latex');
        set(0,'defaultLegendInterpreter','latex');
        set(0,'defaultAxesTickLabelInterpreter','latex');
    end

    set(0,'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',baseFont, ...
          'defaultLineLineWidth',lineLW, ...
          'DefaultAxesBox','off');

    % Color palette (match your Path B feel)
    C.black = [0.10 0.10 0.10];
    C.blue  = [0.00 0.45 0.74];
    C.red   = [0.80 0.25 0.25];
    C.gray  = [0.55 0.55 0.60];
    C.green = [0.20 0.62 0.20];
    C.purp  = [0.42 0.24 0.60];
    C.cyan  = [0.30 0.75 0.80];

    % ======================================================================
    % Fig C1: Input acceleration + basic distribution
    % ======================================================================
    fig1 = figure('Name','Path C — Input Observation', 'Color','w','Units','pixels');
    fig1.Position(3:4) = [980 640];
    tlo1 = tiledlayout(fig1,2,2,'TileSpacing','compact','Padding','compact');

    % Time series (raw)
    ax1 = nexttile(tlo1,1,[1 2]); hold(ax1,'on');
    plot(ax1, t, O, 'Color', C.black, 'LineWidth',1.4, 'DisplayName','a_s(t)');
    applyTheme(ax1, useMinorGrid);
    xlim(ax1,[t(1) t(end)]);
    xlabel(ax1,'Time t [s]');
    ylabel(ax1,'Sprung accel. a_s(t) [m/s^2]');
    title(ax1,'Observation used by ICA (paper: O(t) = a_s(t))');
    niceLegend(ax1,'best');

    % Histogram
    ax2 = nexttile(tlo1,3); hold(ax2,'on');
    histogram(ax2, O, 60, 'Normalization','pdf', 'FaceAlpha',0.9, 'EdgeColor','none');
    applyTheme(ax2, useMinorGrid);
    xlabel(ax2,'a_s [m/s^2]');
    ylabel(ax2,'PDF');
    title(ax2,'Distribution (non-Gaussianity cue)');

    % Embedded channels preview (if enabled)
    ax3 = nexttile(tlo1,4); hold(ax3,'on');
    if use_delay_embedding && M >= 3
        idx_show = round(linspace(1,M,min(5,M)));
        for ii = 1:numel(idx_show)
            plot(ax3, t_eff, Xobs(idx_show(ii),:), 'LineWidth',1.0, ...
                'DisplayName', sprintf('ch %d', idx_show(ii)));
        end
        title(ax3, sprintf('Delay-embedding preview (M=%d, \\tau=%d samples)', M, embed_delay_samples));
        ylabel(ax3,'Embedded accel. [m/s^2]');
        xlabel(ax3,'Time t [s]');
        niceLegend(ax3,'best');
    else
        text(ax3, 0.5, 0.5, 'Delay embedding OFF (single-channel ICA is scale-only)', ...
            'HorizontalAlignment','center', 'Units','normalized', 'Color', C.gray);
        axis(ax3,'off');
    end
    applyTheme(ax3, useMinorGrid);

    export_pub(fig1, fig_dir, 'C01_Observation', fig_format, fig_dpi, save_plots);

    % ======================================================================
    % Fig C2: Whitening diagnostics (covariance heatmaps)
    % ======================================================================
    fig2 = figure('Name','Path C — Whitening Diagnostics', 'Color','w','Units','pixels');
    fig2.Position(3:4) = [980 640];
    tlo2 = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');

    % Cov of centered X
    ax4 = nexttile(tlo2,1);
    imagesc(ax4, covX);
    axis(ax4,'image');
    colorbar(ax4);
    title(ax4,'Covariance of centered observation: cov(X_c)');
    xlabel(ax4,'Channel'); ylabel(ax4,'Channel');
    set(ax4,'YDir','normal');
    applyTheme(ax4, false);

    % Cov of whitened Z
    covZ = (Z*Z.') / max(N_eff-1,1);
    ax5 = nexttile(tlo2,2);
    imagesc(ax5, covZ);
    axis(ax5,'image');
    colorbar(ax5);
    title(ax5,'Covariance after whitening: cov(Z) \approx I');
    xlabel(ax5,'Channel'); ylabel(ax5,'Channel');
    set(ax5,'YDir','normal');
    applyTheme(ax5, false);

    export_pub(fig2, fig_dir, 'C02_Whitening', fig_format, fig_dpi, save_plots);

    % ======================================================================
    % Fig C3: Extracted component in time + spectrum + histogram
    % ======================================================================
    fig3 = figure('Name','Path C — ICA Output', 'Color','w','Units','pixels');
    fig3.Position(3:4) = [980 640];
    tlo3 = tiledlayout(fig3,2,2,'TileSpacing','compact','Padding','compact');

    % Time-domain reconstructed road (ICA source)
    ax6 = nexttile(tlo3,1,[1 2]); hold(ax6,'on');
    plot(ax6, t_ica, zr_t_ica, 'Color', C.green, 'LineWidth',1.6, 'DisplayName','\hat z_r(t) (ICA)');
    applyTheme(ax6, useMinorGrid);
    xlim(ax6, [t_ica(1) t_ica(end)]);
    xlabel(ax6,'Time t [s]');
    ylabel(ax6,'Estimated road \hat z_r(t) [m]');
    title(ax6, sprintf('ICA reconstructed source (scale \\alpha = %.3g), kurtosis = %.3g', ...
        scale_alpha, ica_info.kurt(1)));
    niceLegend(ax6,'best');

    % Spectrum (magnitude)
    ax7 = nexttile(tlo3,3); hold(ax7,'on');
    [f1, P1] = one_sided_psd(zr_t_ica, fs);
    plot(ax7, f1, 10*log10(max(P1, eps)), 'Color', C.blue, 'LineWidth',1.5, 'DisplayName','PSD(\hat z_r)');
    if do_bandlimit_source
        xline(ax7, f_hp, ':', 'HP', 'Color', C.gray, 'LineWidth',1.0, 'HandleVisibility','off');
        xline(ax7, f_lp, ':', 'LP', 'Color', C.gray, 'LineWidth',1.0, 'HandleVisibility','off');
    end
    set(ax7,'XScale','log');
    applyTheme(ax7, useMinorGrid);
    xlabel(ax7,'Frequency [Hz]');
    ylabel(ax7,'PSD [dB re 1 m^2/Hz]');
    title(ax7,'Reconstructed road spectrum');
    niceLegend(ax7,'best');

    % Histogram (non-Gaussianity)
    ax8 = nexttile(tlo3,4); hold(ax8,'on');
    histogram(ax8, zr_t_ica, 60, 'Normalization','pdf', 'FaceAlpha',0.9, 'EdgeColor','none');
    applyTheme(ax8, useMinorGrid);
    xlabel(ax8,'\hat z_r [m]');
    ylabel(ax8,'PDF');
    title(ax8,'ICA source distribution');

    export_pub(fig3, fig_dir, 'C03_ICA_Output', fig_format, fig_dpi, save_plots);

    % ======================================================================
    % Fig C4: Spatial comparison vs ground truth (if available)
    % ======================================================================
    if isfield(cfg,'ground') && isfield(cfg.ground,'x') && isfield(cfg.ground,'zr')
        xg = cfg.ground.x(:);
        zg = cfg.ground.zr(:);

        [x_common, zr_est_resampled, zg_resampled] = ...
            align_on_common_grid(sp.x, sp.zr_x_filt, xg, zg);

        err = zr_est_resampled - zg_resampled;

        rmse = sqrt(mean(err.^2,'omitnan'));
        mae  = mean(abs(err),'omitnan');

        vp = isfinite(zr_est_resampled) & isfinite(zg_resampled);
        if nnz(vp) >= 2 && std(zr_est_resampled(vp)) > 0 && std(zg_resampled(vp)) > 0
            rho = corr(zr_est_resampled(vp), zg_resampled(vp));
        else
            rho = NaN;
        end

        emax = robust_max(abs(err));
        emax = max(emax, eps);
        elims = 1.05*[-emax, emax];

        fig4 = figure('Name','Path C — Spatial Comparison', 'Color','w','Units','pixels');
        fig4.Position(3:4) = [980 640];
        tlo4 = tiledlayout(fig4,2,1,'TileSpacing','compact','Padding','compact');

        ax9 = nexttile(tlo4,1); hold(ax9,'on');
        plot(ax9, xg, zg, '-', 'Color', C.gray, 'LineWidth',1.5, 'DisplayName','Ground truth');
        plot(ax9, sp.x, sp.zr_x_filt, '-', 'Color', C.green, 'LineWidth',1.8, 'DisplayName','ICA recon');
        applyTheme(ax9, useMinorGrid);
        ylabel(ax9,'Elevation z_r(x) [m]');
        title(ax9, sprintf('Spatial profile (RMSE=%.3g m, MAE=%.3g m, \\rho=%.3f)', rmse, mae, rho));
        niceLegend(ax9,'best');

        % --- Bin-wise RMSE over distance (e.g., every 10 m) ---
        ax10 = nexttile(tlo4,2); hold(ax10,'on');
        
        % Choose x axis for binning
        if rmse_bin_use_robust_x
            xb = x_common(:);
            eb = err(:);
        else
            % fallback: use estimate grid directly (less recommended)
            xb = sp.x(:);
            eb = (sp.zr_x_filt(:) - interp1(xg, zg, sp.x(:), 'linear', 'extrap'));
        end
        
        % Define bins
        xmin = min(xb); xmax = max(xb);
        edges = (floor(xmin/rmse_bin_m)*rmse_bin_m) : rmse_bin_m : (ceil(xmax/rmse_bin_m)*rmse_bin_m);
        if numel(edges) < 2
            edges = [xmin, xmax];
        end
        
        rmse_bin = nan(numel(edges)-1, 1);
        x_mid    = nan(numel(edges)-1, 1);
        
        for iBin = 1:(numel(edges)-1)
            inBin = (xb >= edges(iBin)) & (xb < edges(iBin+1)) & isfinite(eb);
            x_mid(iBin) = 0.5*(edges(iBin) + edges(iBin+1));
        
            if nnz(inBin) >= rmse_bin_min_points
                rmse_bin(iBin) = sqrt(mean(eb(inBin).^2, 'omitnan'));
            end
        end
        
        % Plot
        switch lower(rmse_bin_plot_style)
            case 'stairs'
                stairs(ax10, x_mid, rmse_bin, '-', 'Color', C.red, 'LineWidth',1.6, ...
                    'DisplayName', sprintf('RMSE per %.0f m', rmse_bin_m));
            case 'bar'
                bar(ax10, x_mid, rmse_bin, 1.0, 'FaceAlpha',0.75, ...
                    'DisplayName', sprintf('RMSE per %.0f m', rmse_bin_m));
            case 'line'
                plot(ax10, x_mid, rmse_bin, '-', 'Color', C.red, 'LineWidth',1.6, ...
                    'DisplayName', sprintf('RMSE per %.0f m', rmse_bin_m));
            otherwise
                error('Unknown rmse_bin_plot_style: %s', rmse_bin_plot_style);
        end
        
        applyTheme(ax10, useMinorGrid);
        xlabel(ax10,'Distance x [m]');
        ylabel(ax10,'RMSE [m]');
        title(ax10, sprintf('Local RMSE vs distance (bin = %.0f m)', rmse_bin_m));
        niceLegend(ax10,'best');
        linkaxes([ax9,ax10],'x');

        export_pub(fig4, fig_dir, 'C04_Spatial_Comparison', fig_format, fig_dpi, save_plots);
    end

    % ======================================================================
    % Fig C5: ICA convergence summary (iterations + kurtosis)
    % ======================================================================
    fig5 = figure('Name','Path C — ICA Diagnostics', 'Color','w','Units','pixels');
    fig5.Position(3:4) = [980 420];
    tlo5 = tiledlayout(fig5,1,2,'TileSpacing','compact','Padding','compact');

    ax11 = nexttile(tlo5,1); hold(ax11,'on');
    bar(ax11, 1:nComp, ica_info.iter, 'FaceAlpha',0.9);
    applyTheme(ax11, useMinorGrid);
    xlabel(ax11,'Component #');
    ylabel(ax11,'Iterations');
    title(ax11,'Convergence iterations');

    ax12 = nexttile(tlo5,2); hold(ax12,'on');
    bar(ax12, 1:nComp, ica_info.kurt, 'FaceAlpha',0.9);
    yline(ax12, 0, ':', 'Color', C.gray, 'LineWidth',1.0, 'HandleVisibility','off');
    applyTheme(ax12, useMinorGrid);
    xlabel(ax12,'Component #');
    ylabel(ax12,'Excess kurtosis');
    title(ax12,'Non-Gaussianity (kurtosis)');

    export_pub(fig5, fig_dir, 'C05_ICA_Diagnostics', fig_format, fig_dpi, save_plots);
end

fprintf('=== PATH C: done ===\n\n');

%% ------------------------------ HELPERS -----------------------------------
function [Xemb, t_eff, x_eff] = build_delay_embedding(x, t, x_t, M, tau, cropPolicy)
% build_delay_embedding
% Creates an M-channel delay embedding from a single signal x(t):
%   Xemb(1,:) = x(n)
%   Xemb(2,:) = x(n - tau)
%   ...
%   Xemb(M,:) = x(n - (M-1)tau)
%
% cropPolicy:
%   'valid' -> keep only samples where all delays exist (shorter signal)
%   'same'  -> pad with zeros to keep original length (less recommended)

    x = x(:);
    N = numel(x);
    assert(M >= 1 && tau >= 1, 'Embedding: M>=1 and tau>=1 required.');

    maxLag = (M-1)*tau;

    switch lower(cropPolicy)
        case 'valid'
            n0 = 1 + maxLag;
            idx = (n0:N).';
            N_eff = numel(idx);
            Xemb = zeros(M, N_eff);
            for k = 1:M
                lag = (k-1)*tau;
                Xemb(k,:) = x(idx - lag).';
            end
            t_eff = t(idx);
            x_eff = x_t(idx);

        case 'same'
            Xemb = zeros(M, N);
            for k = 1:M
                lag = (k-1)*tau;
                Xtmp = zeros(N,1);
                Xtmp((1+lag):end) = x(1:end-lag);
                Xemb(k,:) = Xtmp.';
            end
            t_eff = t(:);
            x_eff = x_t(:);

        otherwise
            error('Embedding: unknown cropPolicy: %s', cropPolicy);
    end
end

function [Z, Wwhite, Dewhite, covX] = whiten_channels(X)
% whiten_channels
% Centers assumed already done externally if desired; here we compute covariance
% and apply PCA whitening:
%   covX = (X X^T)/(N-1) = E V E^T
%   Wwhite  = V^{-1/2} E^T
%   Dewhite = E V^{1/2}
%   Z       = Wwhite X

    [M, N] = size(X);
    covX = (X*X.') / max(N-1,1);

    % Symmetrize for numerical safety
    covX = 0.5*(covX + covX.');

    [E, V] = eig(covX);
    v = diag(V);

    % Sort eigenvalues descending
    [v, idx] = sort(v, 'descend');
    E = E(:, idx);

    % Regularize tiny/negative eigenvalues
    v = max(v, eps);

    DinvSqrt = diag(1./sqrt(v));
    DSqrt    = diag(sqrt(v));

    Wwhite  = DinvSqrt * E.';   % [M x M]
    Dewhite = E * DSqrt;        % [M x M]

    Z = Wwhite * X;             % [M x N]
end

function w = deflate_orth(w, Wprev)
% deflate_orth
% Gram–Schmidt orthogonalization of w against previously found separating
% vectors (rows of Wprev).
    if isempty(Wprev)
        return;
    end
    for i = 1:size(Wprev,1)
        wi = Wprev(i,:).';
        w  = w - (w.'*wi)*wi;
    end
end

function y = fft_bandlimit(x, fs, f_hp, f_lp)
% fft_bandlimit
% Zero-phase band-limiting by masking FFT bins.
    x = x(:);
    N = numel(x);

    X = fft(x);
    f = (0:N-1)' * (fs/N);

    mask = true(N,1);

    if ~isempty(f_hp) && isfinite(f_hp) && f_hp > 0
        mask = mask & (f >= f_hp) & (f <= fs - f_hp);
    end
    if ~isempty(f_lp) && isfinite(f_lp) && f_lp > 0
        mask = mask & (f <= f_lp | f >= fs - f_lp);
    end

    X(~mask) = 0;
    y = ifft(X, 'symmetric');
end

function y = butter_bandlimit(x, fs, f_hp, f_lp, order)
% butter_bandlimit
% Zero-phase Butterworth band-pass (or HP/LP) using filtfilt.
    x = x(:);
    nyq = fs/2;

    doHP = ~isempty(f_hp) && isfinite(f_hp) && f_hp > 0;
    doLP = ~isempty(f_lp) && isfinite(f_lp) && f_lp > 0 && f_lp < nyq;

    if doHP && doLP && f_hp < f_lp
        [b,a] = butter(order, [f_hp f_lp]/nyq, 'bandpass');
    elseif doHP
        [b,a] = butter(order, f_hp/nyq, 'high');
    elseif doLP
        [b,a] = butter(order, f_lp/nyq, 'low');
    else
        y = x;
        return;
    end

    y = filtfilt(b,a,x);
end

function [f1, P1] = one_sided_psd(x, fs)
% one_sided_psd
% Simple one-sided periodogram PSD (no windowing; good for quick diagnostics).
    x = x(:);
    N = numel(x);
    x = x - mean(x,'omitnan');

    X = fft(x);
    P2 = (abs(X)/N).^2;           % power
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);  % one-sided
    f1 = fs*(0:floor(N/2))'/N;
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
                exportgraphics(fig, [fp '.png'], 'Resolution', dpi);
            case 'tiff'
                exportgraphics(fig, [fp '.tif'], 'Resolution', dpi);
            case 'pdf'
                exportgraphics(fig, [fp '.pdf'], 'ContentType','vector');
            otherwise
                exportgraphics(fig, [fp '.png'], 'Resolution', dpi);
        end
    catch
        switch lower(fmt)
            case 'png'
                print(fig, [fp '.png'], '-dpng', sprintf('-r%d',dpi));
            case 'tiff'
                print(fig, [fp '.tif'], '-dtiff', sprintf('-r%d',dpi));
            case 'pdf'
                set(fig,'PaperPositionMode','auto');
                print(fig, [fp '.pdf'], '-dpdf','-vector');
            otherwise
                print(fig, [fp '.png'], '-dpng', sprintf('-r%d',dpi));
        end
    end
    try
        savefig(fig, [fp '.fig']);
    catch
    end
end

function applyTheme(ax, useMinor)
% Consistent axis styling for all Path C figures.
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

function mu_opt = robustica_opt_mu_kurtosis_real(w, g, Z, targetSign)
% robustica_opt_mu_kurtosis_real
% RobustICA optimal step size for kurtosis contrast (real-valued case).
% Based on Zarzoso & Comon: maximize |K(w + mu g)| by selecting among roots
% of the quartic p(mu) = 0 (derivative of kurtosis along direction). :contentReference[oaicite:6]{index=6}
%
% targetSign:
%   0  -> maximize absolute kurtosis (paper usage)
%   +1 -> target positive kurtosis
%   -1 -> target negative kurtosis

    % Current and direction outputs
    y  = w.' * Z;     % [1 x N]
    yg = g.' * Z;     % [1 x N]
    N  = numel(y);

    % Variables from RobustICA appendix notation (real case simplifies conj/abs)
    a = y.^2;
    b = yg.^2;
    c = y .* yg;
    d = y .* yg;  % Re(y*conj(yg)) for real signals

    % Moments / expectations (sample means)
    Ea   = mean(a);
    Eb   = mean(b);
    Ec   = mean(c);

    h0 = mean(a.^2) - (Ea^2);
    h1 = 4*mean(a.*d) - 4*(Ea*Ec);
    h2 = 4*mean(d.^2) + 2*mean(abs(a).*abs(b)) - 4*(abs(Ec)^2) - 2*(Ea*Eb);
    h3 = 4*mean(abs(b).*d) - 4*(Eb*Ec);
    h4 = mean(b.^2) - (Eb^2);

    i0 = mean(abs(a));
    i1 = 2*mean(d);
    i2 = mean(abs(b));

    % Quartic polynomial p(mu)=sum_{k=0}^4 ak mu^k (from RobustICA appendix) :contentReference[oaicite:7]{index=7}
    a0 = -2*h0*i1 + h1*i0;
    a1 = -4*h0*i2 - h1*i1 + 2*h2*i0;
    a2 = -3*h1*i2 + 3*h3*i0;
    a3 = -2*h2*i2 + h3*i1 + 4*h4*i0;
    a4 = -h3*i2 + 2*h4*i1;

    % Roots of quartic (candidates)
    r = roots([a4 a3 a2 a1 a0]);

    % Candidate step sizes: real parts of roots + include 0
    muCand = [real(r(:)); 0];

    % Evaluate contrast for each candidate and pick best
    bestVal = -Inf;
    mu_opt  = 0;

    for k = 1:numel(muCand)
        mu = muCand(k);
        w_try = w + mu*g;
        w_try = w_try / max(norm(w_try), eps);

        y_try = w_try.' * Z;
        K = kurtosis_excess_real(y_try);

        if targetSign == 0
            val = abs(K);
        elseif targetSign > 0
            val = K;
        else
            val = -K;
        end

        if isfinite(val) && val > bestVal
            bestVal = val;
            mu_opt  = mu;
        end
    end
end

function K = kurtosis_excess_real(y)
% excess kurtosis for real-valued signal y
    y = y(:) - mean(y,'omitnan');
    m2 = mean(y.^2,'omitnan');
    m4 = mean(y.^4,'omitnan');
    if m2 > 0
        K = m4/(m2^2) - 3;
    else
        K = NaN;
    end
end