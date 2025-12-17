
%% Path A: Time-Domain Inversion & Spatial Mapping
% -------------------------------------------------------------------------
% This script combines:
%   3A_a) Two-Step ODE Inversion (Time Domain)
%   3A_b) Map Time-Domain Profile to Spatial Domain & Apply Filter
%
% It:
%   1) Loads recon_cfg.mat and performs the time-domain inversion:
%      a_z^s(t) -> z_s(t), z_u(t) -> z_r(t)
%   2) Saves recon_time_domain.mat
%   3) Loads recon_time_domain.mat and maps z_r(t) -> z_r(x),
%      applies spatial LP (optional), and saves recon_spatial.mat
% -------------------------------------------------------------------------

close all; clear; clc;

%% 3A_a) Two-Step ODE Inversion (Time Domain)
% -------------------------------------------------------------------------
% PURPOSE:
% Implements “Path A” by inverting the quarter-car model directly in the
% time domain using two sequential ODE solves. From the measured sprung
% acceleration a_s(t), it first solves for the unsprung motion z_u(t),
% then solves for the road displacement z_r(t).
%
% CORE LOGIC (How it works):
% 1) Load 'recon_cfg.mat' with a_s(t), fs, V, and vehicle parameters
%    (ms, mu, ks, cs, kt, optional ct).
% 2) De-bias a_s(t) and double-integrate -> z_ṡ(t), z_s(t). High-pass both
%    to suppress drift from numerical integration.
% 3) Solve the sprung-mass ODE for the unsprung displacement:
%       C_s z_u̇ + K_s z_u = M_s a_s + C_s z_ṡ + K_s z_s
%    (1st-order LTI; integrate to obtain z_u(t)).
% 4) Differentiate z_u(t) smoothly to obtain z_u̇(t), z_ü(t) (Savitzky–Golay
%    by default to limit noise amplification).
% 5) Solve the tyre/unsprung ODE for the road profile:
%       C_t z_ṙ + K_t z_r = M_u z_ü + (C_s+C_t) z_u̇ + (K_s+K_t) z_u − (C_s z_ṡ + K_s z_s)
%    (1st-order LTI; integrate to obtain z_r(t)).
% 6) Save all time-domain signals (sprung, unsprung, z_r) for mapping to the
%    spatial domain in the next pipeline step.
%
% OUTPUTS:
% - '00_Outputs/03_ReconstructionCore/PathA_TimeDomain/recon_time_domain.mat'
% -------------------------------------------------------------------------

%% ========================== HYPERPARAMETERS ==========================
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};

% (2) This script’s output folder (consistent naming)
save_folder = fullfile('00_Outputs', '03_ReconstructionCore', 'PathA_TimeDomain');
fig_dir     = fullfile(save_folder,'figs');

% (3) Displacement drift control (actual filtering performed here)
hp_disp_fc_hz   = [];                   % [] -> use value from cfg; else override (Hz)
hp_disp_order   = [];                   % [] -> use value from cfg; else override

% (4) Smoothing settings for derivatives (stabilize z_u̇, z_ü)
use_sgolay      = true;                 % prefer Savitzky–Golay if available
sgolay_win_sec  = 0.5;                  % window length in seconds (auto-rounded to odd samples)
sgolay_poly     = 3;                    % polynomial order (2–5 typical)

% (5) Optional temporal pre-filtering of z_r(t) before spatial mapping
apply_temporal_prefilter = false;       % usually OFF (main LP is spatial)
lp_zr_fc_hz     = 7.0;
lp_zr_order     = 4;

% (6) Plots & saving
make_plots      = true;
save_plots      = true;
fig_dpi         = 300;

%% ----------------------------- Load config -----------------------------
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file') == 2
        cfg_path = cfg_candidates{k};
        break;
    end
end
assert(~isempty(cfg_path), 'Config file not found in expected locations.');

S = load(cfg_path);
assert(isfield(S,'cfg') && isstruct(S.cfg), 'recon_cfg.mat must contain struct cfg');
cfg = S.cfg;

% Ensure output folders exist
if ~exist(save_folder,'dir')
    mkdir(save_folder);
end
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Pull essentials (robust to missing ct/x_t)
t   = cfg.t(:);
fs  = cfg.fs;
dt  = cfg.dt;
V   = cfg.V;
ms  = cfg.ms; mu = cfg.mu; ks = cfg.ks; cs = cfg.cs; kt = cfg.kt;
ct  = 0;
if isfield(cfg,'ct') && ~isempty(cfg.ct)
    ct = cfg.ct;
end
x_t = [];
if isfield(cfg,'x_t') && ~isempty(cfg.x_t)
    x_t = cfg.x_t(:);
end

% Guards / assertions
N = numel(t);
assert(N>=3,'Time vector t is too short.');
assert(all(isfinite([ms mu ks cs kt])),'Vehicle parameters contain NaNs/Inf.');

% Guard tiny/zero dampings for ODE stability
cs_eff = max(cs, 1e-9);
ct_eff = max(ct, 1e-9);

% HP drift control from cfg unless overridden
if isempty(hp_disp_fc_hz)
    hp_disp_fc_hz = getfield_with_default(cfg,'hp_disp_fc_hz', 0.2);
end
if isempty(hp_disp_order)
    hp_disp_order = getfield_with_default(cfg,'hp_disp_order', 2);
end

% Signals
assert(isfield(cfg,'az_s'), 'cfg.az_s (sprung acceleration) missing.');
azs = cfg.az_s(:);                         % sprung vertical acceleration [m/s^2]
assert(numel(azs) == N, 'az_s length must match t.');

% Resonances (for reference/labels only)
f_wh = (1/(2*pi)) * sqrt(max(kt,0)/max(mu,eps));   % wheel-hop (unsprung)
f_sp = (1/(2*pi)) * sqrt(max(ks,0)/max(ms,eps));   % body

%% ----------------------- Step 1: Preprocess a_z^s(t) -----------------------
% De-bias acceleration
azs = azs - mean(azs,'omitnan');

% Double integration (raw) using cumulative trapezoid
vzs_raw = cumtrapz(t, azs);           % z_ṡ (raw)
zss_raw = cumtrapz(t, vzs_raw);       % z_s (raw, drifty)

% High-pass displacement & velocity to control integration drift
[zss, vzs] = hp_disp(zss_raw, vzs_raw, fs, hp_disp_fc_hz, hp_disp_order);

%% ----------------------- Step 2: Solve for z_u(t) --------------------------
% C_s*z_u̇ + K_s*z_u = M_s*a_z^s + C_s*z_ṡ + K_s*z_s  -> first-order LTI ODE
% State: z_u
% dz_u/dt = -(K_s/C_s) * z_u + (1/C_s) * f_s
fs_in = ms*azs + cs*vzs + ks*zss;   % RHS [N]

A1 = -(ks/cs_eff);
B1 =  (1/cs_eff);
C1 =  1;
D1 =  0;

sys1 = ss(A1, B1, C1, D1);          % continuous-time system
zu   = lsim(sys1, fs_in, t);        % unsprung displacement [m]

% Smooth derivative estimates (z_u̇, z_ü)
vzu  = smooth_derivative(zu, t, use_sgolay, sgolay_win_sec, sgolay_poly);
azu  = smooth_derivative(vzu, t, use_sgolay, sgolay_win_sec, sgolay_poly);

%% ----------------------- Step 3: Solve for z_r(t) --------------------------
% C_t*z_ṙ + K_t*z_r = M_u*z_ü + (C_s+C_t)z_u̇ + (K_s+K_t)z_u - (C_s*z_ṡ + K_s*z_s)
rhs = mu*azu + (cs+ct_eff)*vzu + (ks+kt)*zu - (cs*vzs + ks*zss);   % [N]

A2 = -(kt/ct_eff);
B2 =  (1/ct_eff);
C2 =  1;
D2 =  0;

sys2 = ss(A2, B2, C2, D2);
zr_t = lsim(sys2, rhs, t);           % reconstructed road profile (time) [m]

% Optional temporal prefilter (usually off; main low-pass is spatial later)
if apply_temporal_prefilter
    zr_t = lp_temporal(zr_t, fs, lp_zr_fc_hz, lp_zr_order);
end
zr_dot = smooth_derivative(zr_t, t, use_sgolay, sgolay_win_sec, sgolay_poly);

%% ----------------------- QA: Equation balance checks -----------------------
% Check sprung-mass balance used to generate f_s (sanity)
lhs_sprung = ms*azs + cs*vzs + ks*zss;
rhs_sprung = fs_in;

% Check tyre/unsprung balance (C_t*z_ṙ + K_t*z_r) ≟ RHS
lhs_tyre = ct_eff*zr_dot + kt*zr_t;
rhs_tyre = rhs;
residual_tyre = lhs_tyre - rhs_tyre;

%% --------------------------------- Plots -----------------------------------
if make_plots
    % Global styling
    baseFont = 12; axesLW = 1.0; lineLW = 1.4;

    % Color palette (consistent across figures)
    C.az      = [0.00 0.45 0.74];
    C.v_raw   = [0.72 0.82 1.00];
    C.v_hp    = [0.00 0.45 0.74];
    C.z_raw   = [0.75 0.92 0.75];
    C.z_hp    = [0.20 0.65 0.30];
    C.zu      = [0.70 0.00 0.80];
    C.vzu     = [0.85 0.33 0.10];
    C.azu     = [0.10 0.10 0.10];
    C.lhs     = [0.00 0.45 0.74];
    C.rhs     = [0.85 0.33 0.10];
    C.resid   = [0.10 0.10 0.10];

    % =========================================================================
    % Figure 1 — Time-domain pipeline
    % =========================================================================
    f1 = mkfig('Two-Step ODE Inversion — time pipeline', [100 80 1100 760]);
    tlo1 = tiledlayout(f1,4,1,'TileSpacing','compact','Padding','compact');
    title(tlo1,'Time-domain reconstruction pipeline','FontWeight','bold');

    % 1) Sprung acceleration
    ax1 = nexttile(tlo1,1);
    plot(ax1, t, azs, 'Color', C.az, 'LineWidth', lineLW);
    yline(ax1,0,':','Color',[0.5 0.5 0.5]);
    grid(ax1,'on');
    xlabel(ax1,'Time (s)');
    ylabel(ax1,'m/s^2');
    title(ax1,'Sprung acceleration   a_z^s(t)');
    formatAxes(baseFont,axesLW);

    % 2) Integrated velocity
    ax2 = nexttile(tlo1,2); hold(ax2,'on');
    plot(ax2, t, vzs_raw, 'Color', C.v_raw, 'LineWidth', 0.9, ...
         'DisplayName','raw');
    plot(ax2, t, vzs,     'Color', C.v_hp,  'LineWidth', lineLW, ...
         'DisplayName','HP-filtered');
    grid(ax2,'on');
    xlabel(ax2,'Time (s)');
    ylabel(ax2,'m/s');
    title(ax2,'Integrated body velocity   z_ṡ(t)');
    legend(ax2,'Location','best','Box','off');
    formatAxes(baseFont,axesLW);

    % 3) Integrated displacement
    ax3 = nexttile(tlo1,3); hold(ax3,'on');
    plot(ax3, t, zss_raw, 'Color', C.z_raw, 'LineWidth', 0.9, ...
         'DisplayName','raw');
    plot(ax3, t, zss,     'Color', C.z_hp,  'LineWidth', lineLW, ...
         'DisplayName','HP-filtered');
    grid(ax3,'on');
    xlabel(ax3,'Time (s)');
    ylabel(ax3,'m');
    title(ax3,'Integrated body displacement   z_s(t)');
    legend(ax3,'Location','best','Box','off');
    formatAxes(baseFont,axesLW);

    % 4) Unsprung motion
    ax4 = nexttile(tlo1,4); hold(ax4,'on');
    plot(ax4, t, zu,  'Color', C.zu,  'LineWidth', lineLW, ...
         'DisplayName','z_u (m)');
    plot(ax4, t, vzu, 'Color', C.vzu, 'LineWidth', 1.1, ...
         'DisplayName','z_u̇ (m/s)');
    plot(ax4, t, azu, 'Color', C.azu, 'LineWidth', 1.1, ...
         'DisplayName','z_ü (m/s^2)');
    grid(ax4,'on');
    xlabel(ax4,'Time (s)');
    ylabel(ax4,'SI units');
    title(ax4,'Unsprung response   z_u(t), z_u̇(t), z_ü(t)');
    legend(ax4,'Location','best','Box','off');
    formatAxes(baseFont,axesLW);

    linkaxes([ax1 ax2 ax3 ax4],'x');
    maybe_export(f1, fullfile(fig_dir,'A01_time_pipeline.png'), fig_dpi, save_plots);

    % =========================================================================
    % Figure 2 — ODE balance checks
    % =========================================================================
    f2 = mkfig('Two-Step ODE Inversion — equation balance', [130 120 1100 760]);
    tlo2 = tiledlayout(f2,3,1,'TileSpacing','compact','Padding','compact');
    title(tlo2,'Equation balance diagnostics','FontWeight','bold');

    % 1) Sprung equation
    ax1b = nexttile(tlo2,1); hold(ax1b,'on');
    plot(ax1b, t, lhs_sprung, 'Color', C.lhs, 'LineWidth', lineLW, ...
        'DisplayName','LHS');
    plot(ax1b, t, rhs_sprung, '--', 'Color', C.rhs, 'LineWidth', lineLW, ...
        'DisplayName','RHS (=f_s)');
    grid(ax1b,'on');
    xlabel(ax1b,'Time (s)');
    ylabel(ax1b,'N');
    title(ax1b,'Sprung equation:  M_s a_z^s + C_s z_ṡ + K_s z_s');
    legend(ax1b,'Location','best','Box','off');
    formatAxes(baseFont,axesLW);

    % 2) Tyre / unsprung equation
    ax2b = nexttile(tlo2,2); hold(ax2b,'on');
    plot(ax2b, t, lhs_tyre, 'Color', C.lhs, 'LineWidth', lineLW, ...
        'DisplayName','LHS');
    plot(ax2b, t, rhs_tyre, '--', 'Color', C.rhs, 'LineWidth', lineLW, ...
        'DisplayName','RHS');
    grid(ax2b,'on');
    xlabel(ax2b,'Time (s)');
    ylabel(ax2b,'N');
    title(ax2b,'Tyre/unsprung equation:  C_t z_ṙ + K_t z_r');
    legend(ax2b,'Location','best','Box','off');
    formatAxes(baseFont,axesLW);

    % 3) Tyre residual
    ax3b = nexttile(tlo2,3); hold(ax3b,'on');
    plot(ax3b, t, residual_tyre, 'Color', C.resid, 'LineWidth', 1.0);
    yline(ax3b,0,':','Color',[0.4 0.4 0.4]);
    grid(ax3b,'on');
    xlabel(ax3b,'Time (s)');
    ylabel(ax3b,'N');
    title(ax3b,'Tyre equation residual:  (C_t z_ṙ + K_t z_r) - RHS');
    formatAxes(baseFont,axesLW);

    linkaxes([ax1b ax2b ax3b],'x');
    maybe_export(f2, fullfile(fig_dir,'A02_balance_checks.png'), fig_dpi, save_plots);

    % =========================================================================
    % Figure 3 — Spectral diagnostics
    % =========================================================================
    f3 = mkfig('Two-Step ODE Inversion — spectral diagnostics', [160 140 1100 740]);
    tlo3 = tiledlayout(f3,3,1,'TileSpacing','compact','Padding','compact');
    title(tlo3,'Frequency-domain diagnostics (Welch PSD)','FontWeight','bold');

    % 1) PSD of a_z^s
    ax1c = nexttile(tlo3,1);
    [Paz, Faz] = welch_psd(azs, fs);
    if ~isempty(Paz)
        semilogy(ax1c, Faz, Paz, 'Color', C.az, 'LineWidth', 1.2);
        grid(ax1c,'on');
        xlabel(ax1c,'Frequency (Hz)');
        ylabel(ax1c,'PSD (m^2/s^4/Hz)');
        title(ax1c,'PSD of sprung acceleration  a_z^s(t)');
        add_vline(f_sp, '--', 'body f_s_p');
        add_vline(f_wh, '--', 'wheel-hop f_w_h');
        formatAxes(baseFont,axesLW);
    else
        text(0.1,0.5,'Signal too short for PSD','Units','normalized'); axis(ax1c,'off');
    end

    % 2) PSD of z_u
    ax2c = nexttile(tlo3,2);
    [Pzu, Fzu] = welch_psd(zu, fs);
    if ~isempty(Pzu)
        semilogy(ax2c, Fzu, Pzu, 'Color', C.zu, 'LineWidth', 1.2);
        grid(ax2c,'on');
        xlabel(ax2c,'Frequency (Hz)');
        ylabel(ax2c,'PSD (m^2/Hz)');
        title(ax2c,'PSD of unsprung displacement  z_u(t)');
        add_vline(f_sp, '--', 'body f_s_p');
        add_vline(f_wh, '--', 'wheel-hop f_w_h');
        formatAxes(baseFont,axesLW);
    else
        text(0.1,0.5,'Signal too short for PSD','Units','normalized'); axis(ax2c,'off');
    end

    % 3) PSD of reconstructed z_r(t)
    ax3c = nexttile(tlo3,3);
    [Pzr, Fzr] = welch_psd(zr_t, fs);
    if ~isempty(Pzr)
        semilogy(ax3c, Fzr, Pzr, 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
        grid(ax3c,'on');
        xlabel(ax3c,'Frequency (Hz)');
        ylabel(ax3c,'PSD (m^2/Hz)');
        title(ax3c,'PSD of reconstructed road   z_r(t)  (temporal)');
        add_vline(f_wh/3, '-', 'suggested LP (f_w_h/3)');
        formatAxes(baseFont,axesLW);
    else
        text(0.1,0.5,'Signal too short for PSD','Units','normalized'); axis(ax3c,'off');
    end

    maybe_export(f3, fullfile(fig_dir,'A03_spectral_diagnostics.png'), fig_dpi, save_plots);

    % =========================================================================
    % Figure 4 — Optional comparison vs ground truth (if available)
    % =========================================================================
    if isfield(cfg,'ground') && isfield(cfg.ground,'x') && isfield(cfg.ground,'zr')
        if isempty(x_t), x_t = (t - t(1))*V; end
        zr_truth_t = interp1(cfg.ground.x(:), cfg.ground.zr(:), x_t, 'linear', 'extrap');

        f4 = mkfig('Two-Step ODE Inversion — z_r(t) vs ground (if available)', ...
                   [160 140 1100 520]);
        tlo4 = tiledlayout(f4,2,1,'TileSpacing','compact','Padding','compact');
        title(tlo4,'Reconstructed vs. ground-truth road profile','FontWeight','bold');

        % 1) Profiles
        ax1d = nexttile(tlo4,1); hold(ax1d,'on');
        dist = (t - t(1))*V;
        plot(ax1d, dist, zr_truth_t, 'k', 'LineWidth', lineLW, ...
            'DisplayName','Ground truth (spatial→time)');
        plot(ax1d, dist, zr_t,       'Color', C.az, 'LineWidth', 1.1, ...
            'DisplayName','Reconstructed z_r(t)');
        grid(ax1d,'on');
        xlabel(ax1d,'Distance (m)');
        ylabel(ax1d,'z_r (m)');
        legend(ax1d,'Location','best','Box','off');
        title(ax1d,'Temporal profile vs mapped spatial truth');
        formatAxes(baseFont,axesLW);

        % 2) Error
        ax2d = nexttile(tlo4,2); hold(ax2d,'on');
        err = zr_t - zr_truth_t;
        plot(ax2d, dist, err, 'Color', [0.80 0.00 0.00], 'LineWidth', 1.1);
        yline(ax2d,0,':','Color',[0.4 0.4 0.4]);
        grid(ax2d,'on');
        xlabel(ax2d,'Distance (m)');
        ylabel(ax2d,'Error (m)');
        rmse = sqrt(mean(err.^2,'omitnan'));
        mae  = mean(abs(err),'omitnan');
        title(ax2d, sprintf('Error (RMSE = %.3g m, MAE = %.3g m)', rmse, mae));
        formatAxes(baseFont,axesLW);

        linkaxes([ax1d ax2d],'x');
        maybe_export(f4, fullfile(fig_dir,'A04_temporal_vs_truth.png'), fig_dpi, save_plots);
    end
end

%% ----------------------------- Save time-domain outputs -------------------
recon              = struct();
recon.t            = t;
recon.zr_t         = zr_t(:);
recon.zr_dot       = zr_dot(:);
recon.zu           = zu(:);
recon.vzu          = vzu(:);
recon.azu          = azu(:);
recon.zss          = zss(:);
recon.vzs          = vzs(:);
recon.azs          = azs(:);
recon.meta         = struct();
recon.meta.fs      = fs;
recon.meta.dt      = dt;
recon.meta.V       = V;
recon.meta.recon_method = 'Path A: Two-Step ODE Inversion (Time Domain)';
recon.meta.hp_disp_fc_hz = hp_disp_fc_hz;
recon.meta.hp_disp_order = hp_disp_order;
recon.meta.deriv_sgolay = struct('used',use_sgolay,'win_sec',sgolay_win_sec,'poly',sgolay_poly);
recon.meta.temporal_prefilter = struct('applied',apply_temporal_prefilter,'fc_hz',lp_zr_fc_hz,'order',lp_zr_order);
if isfield(cfg,'dx')
    recon.meta.dx = cfg.dx;
end
if isfield(cfg,'fs_spatial')
    recon.meta.fs_spatial = cfg.fs_spatial;
end
if isfield(cfg,'nyq_spatial')
    recon.meta.nyq_spatial = cfg.nyq_spatial;
end
if isfield(cfg,'param_source')
    recon.meta.param_source = cfg.param_source;
end
recon.meta.paths.cfg_path = cfg_path;

out_primary = fullfile(save_folder,'recon_time_domain.mat');
save(out_primary,'recon','-v7.3');
fprintf('Saved: %s\n', out_primary);

%% 3A_b) Map Time-Domain Profile to Spatial Domain & Apply Filter
% -------------------------------------------------------------------------
% PURPOSE:
% Converts time-domain reconstructed road profile z_r(t) into spatial
% domain z_r(x), applies spatial low-pass filter, and saves results.
%
% INPUTS:
% - recon_cfg.mat
% - recon_time_domain.mat (just written above)
%
% OUTPUT:
% - recon_spatial.mat in 00_Outputs/03_ReconstructionCore/SpatialMapping
% -------------------------------------------------------------------------

%% ========================== HYPERPARAMETERS ==========================
% (1) Input discovery
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};

recon_candidates = { ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'PathA_TimeDomain', 'recon_time_domain.mat'), ...
};

% (2) This script’s output folder (your per-script organization)
save_folder_spatial = fullfile('00_Outputs', '03_ReconstructionCore', 'SpatialMapping');
fig_dir_spatial     = fullfile(save_folder_spatial,'figs');

% (3) Spatial grid policy
%     'ground_truth' -> use cfg.ground.x (recommended for later comparisons)
%     'uniform'      -> build uniform grid with dx = V/fs
spatial_grid_policy = 'uniform';   % 'ground_truth' | 'uniform'

% (4) Spatial low-pass (paper default: one-third of wheel-hop, set in cfg.nc_spatial)
%     Leave [] to use cfg values; override here to try alternatives.
lp_nc_cycles_per_m = [];    % [] -> use cfg.nc_spatial; else scalar in cycles/m
lp_order           = [];    % [] -> use cfg.lp_rec_order; typical 4

% (4.5) Toggle spatial filtering
apply_spatial_filter = false;     % true = filter, false = no filtering

% (5) (Optional) DC handling for profile after filtering
restore_mean_after_filter = false;  % add back raw mean after LP (keeps alignment)

% (6) Plots & saving
make_plots_spatial = true;
export_png = true;
% =====================================================================

%% -------------------------- Locate & load inputs --------------------------
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file') == 2
        cfg_path = cfg_candidates{k};
        break;
    end
end
assert(~isempty(cfg_path), 'recon_cfg.mat not found in expected locations.');

recon_path = '';
for k = 1:numel(recon_candidates)
    if exist(recon_candidates{k},'file') == 2
        recon_path = recon_candidates{k};
        break;
    end
end
assert(~isempty(recon_path), 'recon_time_domain.mat not found in expected locations.');

C = load(cfg_path);
cfg = C.cfg;
RT = load(recon_path);
recon = RT.recon;

% Ensure output folders
if ~exist(save_folder_spatial,'dir')
    mkdir(save_folder_spatial);
end
if ~exist(fig_dir_spatial,'dir')
    mkdir(fig_dir_spatial);
end

%% ------------------------------- Essentials -------------------------------
t      = recon.t(:);          % s
zr_t   = recon.zr_t(:);       % m
fs     = cfg.fs;              % Hz
V      = cfg.V;               % m/s
x_time = cfg.x_t(:);          % m (time->distance mapping from init script)

% Choose spatial grid
switch lower(spatial_grid_policy)
    case 'ground_truth'
        x_grid = cfg.ground.x(:);
        grid_label = 'ground_truth';
    case 'uniform'
        x_grid = (t - t(1)) * V;    % same support as time signal
        grid_label = 'uniform_dx = V/fs';
    otherwise
        error('Unknown spatial_grid_policy="%s"', spatial_grid_policy);
end

% Interpolate zr_t onto the chosen spatial grid
zr_x_raw = interp1(x_time, zr_t, x_grid, 'linear', 'extrap');
zr_x_raw(isnan(zr_x_raw)) = 0;

% Spatial sampling info
dx         = mean(diff(x_grid),'omitnan');
fs_spatial = 1 / max(dx, eps);       % samples per meter
nyq_sp     = fs_spatial / 2;         % cycles per meter

% Spatial LP cutoff (cycles/m) and order
nc  = ifelse_empty(lp_nc_cycles_per_m, cfg.nc_spatial);
ord = ifelse_empty(lp_order,           cfg.lp_rec_order);

% Spatial low-pass (optional)
if apply_spatial_filter
    Wn = clamp01(nc / nyq_sp);                 % normalized cutoff
    [b_lp, a_lp] = butter(max(1,ord), Wn, 'low');
    zr_x_filt = filtfilt(b_lp, a_lp, zr_x_raw);

    % (optional) restore mean to match raw
    if restore_mean_after_filter
        mraw = mean(zr_x_raw,'omitnan');
        mf = mean(zr_x_filt,'omitnan');
        zr_x_filt = zr_x_filt + (mraw - mf);
    end
else
    zr_x_filt = zr_x_raw;  % pass-through
end

%% --------------------------------- Plots ----------------------------------
if make_plots_spatial
    % ---------- Styling ----------
    baseFont = 12; axesLW = 0.9; lineLW = 1.6;
    set(0,'defaultAxesFontName','Calibri','defaultAxesFontSize',baseFont);
    set(0,'defaultLineLineWidth',lineLW);
    set(0,'DefaultAxesBox','on');

    % Color palette
    Cc.gray  = [0.55 0.55 0.55];
    Cc.blue  = [0.00 0.45 0.74];
    Cc.light = [0.80 0.86 1.00];
    Cc.red   = [0.85 0.33 0.10];
    Cc.dk    = [0.15 0.15 0.15];

    % ---------------------------------------------------------------------
    % Fig 1 — Mapping overview: time → spatial
    % ---------------------------------------------------------------------
    f1s = mkfig('Spatial mapping overview', [120 90 1050 420]);
    tl1 = tiledlayout(f1s,2,1,'TileSpacing','compact','Padding','compact');

    % Time domain
    ax = nexttile(tl1,1); hold(ax,'on');
    plot(ax, t, zr_t, 'Color', Cc.dk);
    grid(ax,'on'); xlabel(ax,'Time (s)'); ylabel(ax,'Elevation (m)');
    title(ax,'Reconstructed road — time domain z_r(t)');

    % Spatial domain (raw)
    ax = nexttile(tl1,2); hold(ax,'on');
    plot(ax, x_grid, zr_x_raw, 'Color', Cc.blue);
    grid(ax,'on'); xlabel(ax,'Distance (m)'); ylabel(ax,'Elevation (m)');
    title(ax, sprintf('Reconstructed road — spatial domain z_r(x)  [%s grid]   (dx=%.3g m, f_s=%.3g samp/m)', ...
        strrep(grid_label,'_','\_'), dx, fs_spatial));
    ylim(ax,'auto');

    if export_png
        exportgraphics(f1s, fullfile(fig_dir_spatial,'07_mapping_overview.png'), 'Resolution', 200);
    end

    % ---------------------------------------------------------------------
    % Fig 2 — Full profile + two auto-zoom panels (windows highlighted)
    % ---------------------------------------------------------------------
    f2s = mkfig('Spatial profile: raw vs filtered', [140 110 1150 640]);
    tl2 = tiledlayout(f2s,2,2,'TileSpacing','compact','Padding','compact');

    % Top: full length, draw zoom boxes
    ax1s = nexttile(tl2,[1 2]); hold(ax1s,'on');
    if apply_spatial_filter
        plot(ax1s, x_grid, zr_x_raw,  'Color', Cc.light, 'LineWidth', 0.9, 'DisplayName','raw');
        plot(ax1s, x_grid, zr_x_filt, 'Color', Cc.blue,  'LineWidth', 1.5, 'DisplayName','filtered');
        title(ax1s, sprintf('z_r(x) — full length   (LP: n_c = %.3f cyc/m, order %d)', nc, ord));
        legend(ax1s,'Location','best');
    else
        plot(ax1s, x_grid, zr_x_raw,  'Color', Cc.blue,  'LineWidth', 1.5, 'DisplayName','raw');
        title(ax1s, 'z_r(x) — full length (spatial LP disabled)');
        legend(ax1s,'Location','best');
    end
    grid(ax1s,'on'); xlabel(ax1s,'Distance (m)'); ylabel(ax1s,'Elevation (m)');
    ylim(ax1s,'auto');

    % Choose two zoom centers from the filtered (or raw) profile
    sig_for_zoom = zr_x_filt;
    if ~apply_spatial_filter, sig_for_zoom = zr_x_raw; end
    [idx1, idx2] = choose_zoom_centers(sig_for_zoom, fs_spatial);

    halfW_m = 20;                         % half window width [m]
    halfW   = max(1, round(halfW_m*fs_spatial));

    % Zoom 1
    idx = max(1,idx1-halfW):min(numel(x_grid),idx1+halfW);
    ax2s = nexttile(tl2); hold(ax2s,'on');
    if apply_spatial_filter
        plot(ax2s, x_grid(idx), zr_x_raw(idx),  'Color', Cc.light, 'LineWidth', 1.0, 'DisplayName','raw');
        plot(ax2s, x_grid(idx), zr_x_filt(idx), 'Color', Cc.blue,  'LineWidth', 1.6, 'DisplayName','filtered');
    else
        plot(ax2s, x_grid(idx), zr_x_raw(idx),  'Color', Cc.blue,  'LineWidth', 1.6, 'DisplayName','raw');
    end
    grid(ax2s,'on'); xlabel(ax2s,'x (m)'); ylabel(ax2s,'z_r (m)');
    title(ax2s, sprintf('Zoom 1 (center @ %.1f m)', x_grid(idx1)));
    legend(ax2s,'Location','best');
    draw_zoom_box(ax1s, x_grid(idx([1 end])), ylim(ax1s), Cc.red, 0.06);

    % Zoom 2
    idx = max(1,idx2-halfW):min(numel(x_grid),idx2+halfW);
    ax3s = nexttile(tl2); hold(ax3s,'on');
    if apply_spatial_filter
        plot(ax3s, x_grid(idx), zr_x_raw(idx),  'Color', Cc.light, 'LineWidth', 1.0, 'DisplayName','raw');
        plot(ax3s, x_grid(idx), zr_x_filt(idx), 'Color', Cc.blue,  'LineWidth', 1.6, 'DisplayName','filtered');
    else
        plot(ax3s, x_grid(idx), zr_x_raw(idx),  'Color', Cc.blue,  'LineWidth', 1.6, 'DisplayName','raw');
    end
    grid(ax3s,'on'); xlabel(ax3s,'x (m)'); ylabel(ax3s,'z_r (m)');
    title(ax3s, sprintf('Zoom 2 (center @ %.1f m)', x_grid(idx2)));
    legend(ax3s,'Location','best');
    draw_zoom_box(ax1s, x_grid(idx([1 end])), ylim(ax1s), Cc.red, 0.06);

    if export_png
        exportgraphics(f2s, fullfile(fig_dir_spatial,'08_profile_full_and_zooms.png'), 'Resolution', 220);
    end

    % ---------------------------------------------------------------------
    % Fig 3 — Spatial spectrum (ASD) with wavelength axis and cutoff
    % ---------------------------------------------------------------------
    [Praw, Fx] = welch_psd_spatial(zr_x_raw, fs_spatial);
    ASD_raw = sqrt(Praw);
    f3s = mkfig('Spatial Spectrum (ASD)', [170 140 1100 520]);
    axS = axes(f3s); hold(axS,'on'); set(axS,'XScale','log','YScale','log','LineWidth',axesLW);
    loglog(axS, Fx, ASD_raw, 'Color', Cc.gray, 'LineWidth', 1.1, 'DisplayName','raw');

    if apply_spatial_filter
        [Pfilt, Fx2] = welch_psd_spatial(zr_x_filt, fs_spatial);
        ASD_filt = sqrt(Pfilt);
        loglog(axS, Fx2, ASD_filt, 'Color', Cc.blue, 'LineWidth', 1.6, 'DisplayName','filtered');
        % Cutoff marker + shaded stopband
        yl = ylim(axS);
        plot(axS, [nc nc], yl, '-', 'Color', Cc.dk, 'LineWidth', 1.0, 'HandleVisibility','off');
        text(nc, yl(2), sprintf('  n_c=%.3f cyc/m', nc), 'Rotation',90, ...
            'VerticalAlignment','top','HorizontalAlignment','left','Color',Cc.dk,'FontSize',9);
        shade_band_x(axS, [nc max(Fx)], yl, [0.9 0.9 0.9], 0.35);
        legend(axS,'Location','southwest');
    else
        legend(axS,'Location','southwest');
    end

    grid(axS,'on'); xlabel(axS,'Spatial frequency f_x (cycles/m)');
    ylabel(axS,'ASD  (m / \surd(cycles/m))');
    title(axS,'Spatial Amplitude Spectral Density (Welch)');

    % Add top wavelength axis: lambda = 1 / f_x
    overlay_wavelength_axis(axS);

    if export_png
        exportgraphics(f3s, fullfile(fig_dir_spatial,'09_spatial_spectrum_asd.png'), 'Resolution', 220);
    end

    % ---------------------------------------------------------------------
    % Fig 4 — Filter magnitude actually applied (only if filtering)
    % ---------------------------------------------------------------------
    if apply_spatial_filter
        [b_lp2, a_lp2] = butter(max(1,ord), clamp01(nc/(fs_spatial/2)), 'low');
        [Hf,w] = freqz(b_lp2, a_lp2, 4096);     % rad/sample
        f_norm = w/(2*pi);                    % cycles/sample
        f_sp   = f_norm * fs_spatial;         % cycles/m

        f4s = mkfig('Spatial LP magnitude', [200 170 900 480]);
        axF = axes(f4s); hold(axF,'on'); set(axF,'XScale','log','LineWidth',axesLW);
        plot(axF, f_sp, 20*log10(abs(Hf)+eps), 'Color', Cc.blue, 'LineWidth', 1.6);
        grid(axF,'on'); xlabel(axF,'Spatial frequency f_x (cycles/m)');
        ylabel(axF,'|H_{LP}(f_x)| (dB)');
        title(axF, sprintf('Applied spatial Butterworth LP — order %d, n_c=%.3f cyc/m (dx=%.4g m, f_{Nyq}=%.3g cyc/m)', ...
            ord, nc, dx, fs_spatial/2));
        xlim(axF,[min(f_sp(f_sp>0)), fs_spatial/2*1.05]);
        yline(axF, -3, ':', 'Color', Cc.gray, 'HandleVisibility','off');
        plot(axF, [nc nc], ylim(axF), '-', 'Color', Cc.dk, 'HandleVisibility','off');
        overlay_wavelength_axis(axF);

        if export_png
            exportgraphics(f4s, fullfile(fig_dir_spatial,'10_filter_magnitude.png'), 'Resolution', 220);
        end
    end
end

%% ----------------------------- Save spatial outputs -----------------------
sp = struct();
sp.x               = x_grid;
sp.zr_x_raw        = zr_x_raw;
sp.zr_x_filt       = zr_x_filt;     % equals raw when filtering disabled
sp.dx              = dx;
sp.fs_spatial      = fs_spatial;
sp.nc_spatial      = nc;
sp.lp_order        = ord;
sp.filter_applied  = logical(apply_spatial_filter);
sp.meta.grid_policy= spatial_grid_policy;
sp.meta.fs         = fs;
sp.meta.V          = V;
sp.meta.cfg_snapshot = cfg;

sp_out = fullfile(save_folder_spatial,'recon_spatial.mat');
save(sp_out,'sp','-v7.3');
fprintf('Saved: %s  (filter_applied = %d)\n', sp_out, sp.filter_applied);


%% =============================== FUNCTIONS ===============================
function val = getfield_with_default(s, name, defaultVal)
    if isfield(s,name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = defaultVal;
    end
end

function [z_hp, vz_hp] = hp_disp(z_raw, vz_raw, fs, fc, order)
    % High-pass the integrated displacement & velocity to remove drift.
    if isempty(fc) || fc<=0 || ~isfinite(fc)
        z_hp = z_raw;
        vz_hp = vz_raw;
        return;
    end
    if fc >= fs/2
        fc = 0.99*(fs/2);
    end
    ord = max(1, round(order));
    Wn = clamp01(fc/(fs/2));
    [b,a] = butter(ord, Wn, 'high');
    z_hp  = filtfilt(b,a, z_raw);
    vz_hp = filtfilt(b,a, vz_raw);
end

function zr_f = lp_temporal(zr, fs, fc, order)
    if isempty(fc) || fc<=0 || ~isfinite(fc)
        zr_f = zr;
        return;
    end
    if fc >= fs/2
        fc = 0.99*(fs/2);
    end
    ord = max(1, round(order));
    Wn = clamp01(fc/(fs/2));
    [b,a] = butter(ord, Wn, 'low');
    zr_f = filtfilt(b,a, zr);
end

function v = smooth_derivative(y, t, useSG, win_sec, poly)
    % Compute dy/dt with smoothing (Savitzky–Golay preferred).
    y = y(:); t = t(:);
    N = numel(y);
    if N < 3 || numel(t) ~= N
        v = zeros(size(y));
        return;
    end
    dt = mean(diff(t),'omitnan');
    if ~isfinite(dt) || dt <= 0
        dt = (t(end)-t(1))/max(N-1,1);
    end

    if useSG && exist('sgolayfilt','file')==2
        % SG window (odd number of samples), clamped to [5, N-1] and odd
        w = max(5, 2*floor((win_sec/dt)/2)+1);
        if mod(w,2) == 0
            w = w + 1;
        end
        w = min(w, N - (mod(N,2) == 0));
        w = max(5, w);
        poly = min(max(round(poly),2),5);
        y_s = sgolayfilt(y, poly, max(w,3));
        v = gradient(y_s, t);
    else
        % Fallback: light moving-average then gradient
        w = max(5, 2*floor((win_sec/dt)/2)+1);
        if mod(w , 2) == 0
            w = w + 1;
        end
        w = min(w, N - (mod(N,2)==0));
        y_s = movmean(y, w, 'Endpoints','shrink');
        v   = gradient(y_s, t);
    end
end

function [Pxx,F] = welch_psd(x, fs)
    % Robust Welch PSD: never use a segment longer than the signal.
    x = x(:) - mean(x(:),'omitnan');
    N = numel(x);
    if N < 8
        Pxx = [];
        F = [];
        return;
    end
    % Choose nperseg ≤ N, about N/4 but at least 128 (if possible)
    nperseg = min(N, max(128, 2^floor(log2(max(32, floor(N/4))))));
    win  = hamming(nperseg);
    nover= floor(nperseg/2);
    nfft = 2^nextpow2(max(nperseg, min(8192, N)));
    [Pxx,F] = pwelch(x, win, nover, nfft, fs, 'onesided');
end

function f = mkfig(name, pos)
    if nargin<2
        pos = [100 100 1000 700];
    end
    f = figure('Name', name, 'Color', 'w', 'Position', pos);
end

function formatAxes(baseFont, axesLW)
    if nargin < 1
        baseFont = 11;
    end
    if nargin < 2
        axesLW = 1.0;
    end
    ax = gca;
    ax.LineWidth = axesLW;
    ax.FontName = 'Calibri';
    ax.FontSize = baseFont;
    grid(ax,'on'); box(ax,'on');
end

function add_vline(x, style, labelStr)
    if isempty(x) || ~isfinite(x)
        return;
    end
    yl = ylim; hold on;
    plot([x x], yl, style, 'Color',[0.2 0.2 0.2], 'LineWidth', 1.0);
    text(x, yl(2), sprintf('  %s (%.2f Hz)', labelStr, x), ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'Rotation',90,'FontSize',9,'Color',[0.2 0.2 0.2]);
end

function y = clamp01(y)
    y = max(min(y,0.999), 1e-6);
end

function maybe_export(figHandle, outPath, dpi, doSave)
    if ~doSave
        return;
    end
    try
        exportgraphics(figHandle, outPath, 'Resolution', dpi);
    catch
        [folder, base, ~] = fileparts(outPath);
        if ~exist(folder,'dir')
            mkdir(folder);
        end
        saveas(figHandle, fullfile(folder, [base '.png']));
    end
end

function [Pxx,F] = welch_psd_spatial(x, fs_spatial)
    % Welch PSD in spatial domain. Units ~ m^2·m (depends on convention).
    x = x(:) - mean(x(:),'omitnan');
    Npts  = numel(x);
    nfft  = max(1024, 2^nextpow2(min(8192, Npts)));
    win   = hamming(round(nfft/2));
    nover = round(numel(win)/2);
    [Pxx, F] = pwelch(x, win, nover, nfft, fs_spatial, 'onesided');
end

function out = ifelse_empty(val, fallback)
    if isempty(val)
        out = fallback;
    else
        out = val;
    end
end

function [xZ, y1Z, y2Z] = middle_window(x, y1, y2, frac)
    % Extract a centered window of length = frac * full length (0<frac<=1).
    if nargin<4 || isempty(frac)
        frac = 0.10;
    end
    frac = min(max(frac, 0.02), 1.0);
    N = numel(x);
    L = max(3, round(frac * N));
    s = floor((N - L)/2) + 1;
    e = s + L - 1;
    xZ  = x(s:e);
    y1Z = y1(s:e);
    y2Z = y2(s:e);
end

function [lo, hi] = robust_ylim(y)
    % 1–99th percentile with a little padding.
    y = y(isfinite(y));
    if isempty(y), lo=-1; hi=1; return; end
    q = prctile(y,[1 99]);
    pad = 0.08*range(q);
    lo = q(1)-pad; hi = q(2)+pad;
end

function [idx1, idx2] = choose_zoom_centers(sig, fs_spatial)
    % Pick two zoom centers using prominent features, well separated.
    s = sig(:) - median(sig(:),'omitnan');
    minSep = max(1, round(10*fs_spatial));           % ≥10 m apart
    prom   = 0.75*std(s,'omitnan'); if ~isfinite(prom) || prom<=0, prom=1e-3; end
    [~,loc,~,pr] = findpeaks(abs(s), 'MinPeakDistance',minSep, 'MinPeakProminence',prom);
    if isempty(loc)
        loc = round(linspace(round(numel(s)*0.25), round(numel(s)*0.75), 2));
        pr  = [1 0];
    end
    [~,ord] = sort(pr,'descend'); loc = loc(ord);
    idx1 = loc(1);
    idx2 = loc(min(2,numel(loc)));
    if idx2==idx1, idx2 = min(numel(s), idx1+minSep); end
end

function draw_zoom_box(ax, xwin, ylims, color, alphaVal)
    % Draw a translucent rectangle indicating the zoom window on a parent axis.
    px = [xwin(1) xwin(2) xwin(2) xwin(1)];
    py = [ylims(1) ylims(1) ylims(2) ylims(2)];
    patch('Parent',ax,'XData',px,'YData',py,'FaceColor',color,'FaceAlpha',alphaVal, ...
          'EdgeColor','none','HandleVisibility','off');
end

function shade_band_x(ax, xlimBand, ylims, color, alphaVal)
    % Shade an x-interval (e.g., stopband) on log or linear axes.
    px = [xlimBand(1) xlimBand(2) xlimBand(2) xlimBand(1)];
    py = [ylims(1) ylims(1) ylims(2) ylims(2)];
    patch('Parent',ax,'XData',px,'YData',py,'FaceColor',color,'FaceAlpha',alphaVal, ...
          'EdgeColor','none','HandleVisibility','off');
end

function overlay_wavelength_axis(axBottom)
    % Add a top x-axis with wavelength λ = 1/f_x (m). Works with log x.
    xl = xlim(axBottom);
    xl(1) = max(xl(1), 1e-3);            % avoid 1/0

    axTop = axes('Position', axBottom.Position, 'Color','none', ...
                 'XAxisLocation','top','YAxisLocation','right','YColor','none', ...
                 'XScale', axBottom.XScale, 'YTick', []);
    xlim(axTop, xl);

    bt = get(axBottom,'XTick');
    bt = bt(bt>0 & bt>=xl(1) & bt<=xl(2));
    if numel(bt) >= 3
        set(axTop,'XTick', bt, 'XTickLabel', arrayfun(@(f)sprintf('%.2g',1./f), bt, 'UniformOutput',false));
    else
        tx = logspace(log10(xl(1)), log10(xl(2)), 6);
        set(axTop,'XTick', tx, 'XTickLabel', arrayfun(@(f)sprintf('%.2g',1./f), tx, 'UniformOutput',false));
    end
    xlabel(axTop,'Wavelength \lambda (m)');
end