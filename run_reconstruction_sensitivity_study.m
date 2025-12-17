
%% 9) Run Parameter Sensitivity Study for the Reconstruction
% -------------------------------------------------------------------------
% PURPOSE:
% This script investigates how sensitive the reconstruction algorithm is to
% inaccuracies in the assumed vehicle parameters. It systematically varies
% each parameter, re-runs the reconstruction, and measures the impact on
% the final error.
%
% CORE LOGIC (How it works):
% 1.  Loads the configuration file ('recon_cfg.mat') to get the input
%     acceleration, ground truth, and baseline vehicle parameters.
%
% 2.  Defines a list of key vehicle parameters to test (e.g., suspension
%     stiffness 'ks', sprung mass 'ms', etc.) and a range of scale factors
%     (e.g., 0.8, 0.9, 1.0, 1.1, 1.2, representing -20%, -10%, etc. error).
%
% 3.  Performs a "one-at-a-time" sensitivity sweep. It loops through each
%     parameter, and for each parameter, it loops through each scale factor.
%
% 4.  In each inner loop, it creates a modified set of vehicle parameters
%     (e.g., ks = baseline_ks * 0.8), and then re-runs the ENTIRE
%     reconstruction pipeline internally (time-domain inversion and spatial
%     mapping).
%
% 5.  It compares the resulting profile to the ground truth and calculates
%     the error (RMSE).
%
% 6.  After all loops are complete, it saves the full set of results,
%     showing how the RMSE changes as each parameter is varied.
%
% 7.  Generates summary plots, including a "tornado plot" that clearly
%     visualizes which parameters have the largest impact on reconstruction
%     accuracy.
%
% INPUTS:
% - '02_Reconstruction Inputs/recon_cfg.mat': The baseline configuration.
%
% OUTPUTS:
% - '09_Sensitivity Study/sensitivity_results.mat': Contains the full
%   dataset from the sensitivity sweep.
% - '09_Sensitivity Study/sensitivity_results.csv': A summary CSV.
% - Figures visualizing the sensitivity results.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ========================== HYPERPARAMETERS ==========================
% (1) Input discovery (new location first, then legacy)
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};

% (2) This script’s output folder
out_root = fullfile('00_Outputs', '04_Validation', 'SensitivityAnalysis');
fig_dir  = fullfile(out_root,'figs');

% (3) Parameters to sweep (toggle any off by setting false)
sweep_flags = struct( ...
    'ks', true, ...   % Suspension stiffness
    'cs', true, ...   % Suspension damping
    'kt', true, ...   % Tyre stiffness
    'mu', true, ...   % Unsprung mass
    'ms', true, ...   % Sprung mass
    'ct', false ...   % Tyre damping (often nominal; enable if you want)
);

% (4) Scale factors for each sweep (multiplicative)
scale_factors = [0.80 0.90 1.00 1.10 1.20];

% (5) Displacement drift control (as used in the paper; can override cfg)
hp_disp_fc_hz_override = [];   % [] => use cfg.hp_disp_fc_hz
hp_disp_order_override = [];   % [] => use cfg.hp_disp_order

% (6) Derivative smoothing (for z_u̇, z_ü)
use_sgolay     = true;
sgolay_win_sec = 0.15;         % seconds
sgolay_poly    = 3;

% (7) Spatial grid & low-pass policy
spatial_grid_policy = 'ground_truth';  % 'ground_truth' | 'uniform'
lp_order_override   = [];              % [] => use cfg.lp_rec_order
restore_mean_after_filter = true;

% (8) Metrics & plots
remove_dc_for_psd = true;
make_plots  = true;
export_png  = true;

%% ------------------------- Load configuration -------------------------
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file')==2, cfg_path = cfg_candidates{k}; break; end
end
assert(~isempty(cfg_path), 'Missing recon_cfg.mat in expected locations.');

C = load(cfg_path); cfg = C.cfg;

if ~exist(out_root,'dir'), mkdir(out_root); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Ground truth (required for sensitivity)
assert(isfield(cfg,'ground') && isfield(cfg.ground,'x') && ~isempty(cfg.ground.x), ...
    'Ground-truth profile not found in cfg.ground.x. Provide ground truth to compare.');

% Pull signals & base params
t   = cfg.t(:);
fs  = cfg.fs;
dt  = cfg.dt;
V   = cfg.V;
x_t = cfg.x_t(:);
azs = cfg.az_s(:);  % sprung acceleration

% Params (baseline)
P0.ms = cfg.ms; P0.mu = cfg.mu; P0.ks = cfg.ks; P0.cs = cfg.cs; P0.kt = cfg.kt; P0.ct = cfg.ct;

% Filter/derivative settings
hp_fc = pickval(hp_disp_fc_hz_override, cfg.hp_disp_fc_hz);
hp_or = pickval(hp_disp_order_override, cfg.hp_disp_order);
lp_or = pickval(lp_order_override,      cfg.lp_rec_order);

% Ground truth grid and profile
x_true  = cfg.ground.x(:);
zr_true = cfg.ground.zr(:);

% Build the evaluation x-grid
switch lower(spatial_grid_policy)
    case 'ground_truth'
        x_grid = x_true;
        grid_label = 'ground_truth';
    case 'uniform'
        x_grid = (t - t(1)) * V;
        grid_label = 'uniform_dx=V/fs';
    otherwise
        error('Unknown spatial_grid_policy="%s"', spatial_grid_policy);
end
dx         = mean(diff(x_grid),'omitnan');
fs_spatial = 1 / max(dx, eps);
nyq_sp     = fs_spatial / 2;

%% ------------------------- Baseline reconstruction -------------------------
base = reconstruct_once(t, azs, fs, V, x_t, x_grid, P0, ...
                        hp_fc, hp_or, use_sgolay, sgolay_win_sec, sgolay_poly, ...
                        lp_or, restore_mean_after_filter);

% Metrics baseline
[Mb, bands_used] = compare_metrics(x_grid, base.zr_x, x_true, zr_true, fs_spatial, remove_dc_for_psd);

%% ------------------------- One-at-a-time sweeps -------------------------
sweep_list = {};
if sweep_flags.ks, sweep_list{end+1} = 'ks'; end
if sweep_flags.cs, sweep_list{end+1} = 'cs'; end
if sweep_flags.kt, sweep_list{end+1} = 'kt'; end
if sweep_flags.mu, sweep_list{end+1} = 'mu'; end
if sweep_flags.ms, sweep_list{end+1} = 'ms'; end
if sweep_flags.ct, sweep_list{end+1} = 'ct'; end

results = struct(); % results.(param).scales, rmse, mae, bias, r, rse_all, etc.

for ip = 1:numel(sweep_list)
    pname = sweep_list{ip};
    P = P0; % reset to baseline each loop
    scales = scale_factors(:)';
    rmse = zeros(size(scales));
    mae  = zeros(size(scales));
    bias = zeros(size(scales));
    rr   = zeros(size(scales));
    rse_all = zeros(size(scales));

    % Store overlay profiles for −20%, 1.0, +20% for later plotting (max three)
    prof_store = struct('scale',[],'x',[],'zr',[]);
    store_idx = [find(abs(scales-0.80)<1e-9,1), find(abs(scales-1.00)<1e-9,1), find(abs(scales-1.20)<1e-9,1)];
    store_idx = store_idx(isfinite(store_idx));

    for is = 1:numel(scales)
        s = scales(is);
        P.(pname) = P0.(pname) * s;

        rec = reconstruct_once(t, azs, fs, V, x_t, x_grid, P, ...
                               hp_fc, hp_or, use_sgolay, sgolay_win_sec, sgolay_poly, ...
                               lp_or, restore_mean_after_filter);

        % Metrics vs truth
        M = compare_metrics(x_grid, rec.zr_x, x_true, zr_true, fs_spatial, remove_dc_for_psd);

        rmse(is) = M.pointwise.rmse_m;
        mae(is)  = M.pointwise.mae_m;
        bias(is) = M.pointwise.bias_m;
        rr(is)   = M.pointwise.pearson_r;
        rse_all(is) = M.spectral.rse_all;

        if any(store_idx == is)
            k = find(store_idx==is);
            prof_store(k).scale = s;
            prof_store(k).x     = x_grid;
            prof_store(k).zr    = rec.zr_x;
        end
    end

    results.(pname).scales  = scales;
    results.(pname).rmse    = rmse;
    results.(pname).mae     = mae;
    results.(pname).bias    = bias;
    results.(pname).r       = rr;
    results.(pname).rse_all = rse_all;
    results.(pname).profiles = prof_store;
end

%% ------------------------- Tornado summary (±20%) -------------------------
% Build ΔRMSE bars relative to baseline for −20% and +20% per parameter.
params_order = sweep_list;
delta_rmse = zeros(numel(params_order), 2);
for i=1:numel(params_order)
    p = params_order{i};
    sc = results.(p).scales;
    [~,i80] = min(abs(sc-0.80));
    [~,i120]= min(abs(sc-1.20));
    delta_rmse(i,1) = results.(p).rmse(i80)  - Mb.pointwise.rmse_m;
    delta_rmse(i,2) = results.(p).rmse(i120) - Mb.pointwise.rmse_m;
end

% Sort by max |ΔRMSE|
[~,ord] = sort(max(abs(delta_rmse),[],2),'descend');
params_sorted = params_order(ord);
delta_sorted  = delta_rmse(ord,:);

%% ------------------------- Save results -------------------------
S = struct();
S.baseline.params = P0;
S.baseline.metrics = Mb;
S.baseline.profile.x = x_grid;
S.baseline.profile.zr = base.zr_x;
S.results = results;
S.tornado.params_sorted = params_sorted;
S.tornado.delta_rmse    = delta_sorted;
S.settings.scale_factors = scale_factors;
S.settings.sweep_flags   = sweep_flags;
S.settings.use_sgolay    = use_sgolay;
S.settings.sgolay_win_sec= sgolay_win_sec;
S.settings.sgolay_poly   = sgolay_poly;
S.settings.hp_fc_hz      = hp_fc;
S.settings.hp_order      = hp_or;
S.settings.lp_order      = lp_or;
S.settings.spatial_grid_policy = spatial_grid_policy;

save(fullfile(out_root,'sensitivity_results.mat'), 'S','-v7.3');
fprintf('Saved: %s\n', fullfile(out_root,'sensitivity_results.mat'));

% CSV summary — one row per (param, scale)
csv_path = fullfile(out_root,'sensitivity_results.csv');
write_sensitivity_csv(csv_path, results);
fprintf('Saved: %s\n', csv_path);

%% ------------------------- Plots -------------------------
if make_plots
    % Fig 1 — RMSE vs scale for each parameter
    f1 = mkfig('Sensitivity — RMSE vs scale', [100 70 1180 760]);
    tlo1 = tiledlayout(f1, numel(params_order), 1, 'TileSpacing','compact','Padding','compact');
    for i=1:numel(params_order)
        p = params_order{i};
        nexttile(tlo1,i);
        plot(results.(p).scales, results.(p).rmse, '-o','LineWidth',1.2,'MarkerSize',4);
        grid on; box on; formatAxes();
        title(sprintf('RMSE vs scale — %s', upper(p)), 'Interpreter','none');
        xlabel('Scale'); ylabel('RMSE (m)');
        yline(S.baseline.metrics.pointwise.rmse_m, 'k--', 'Baseline', 'LabelVerticalAlignment','bottom');
    end
    if export_png, exportgraphics(f1, fullfile(fig_dir,'20_sensitivity_lines.png'), 'Resolution', 180); end

    % Fig 2 — Tornado chart (ΔRMSE at ±20%)
    f2 = mkfig('Sensitivity — Tornado ΔRMSE (±20%)', [140 100 980 620]);
    barh_cool(params_sorted, delta_sorted);
    if export_png, exportgraphics(f2, fullfile(fig_dir,'21_tornado_rmse.png'), 'Resolution', 180); end

    % Fig 3 — Overlay for most sensitive parameter (−20%, 1.0, +20%)
    if ~isempty(params_sorted)
        ptop = params_sorted{1};
        prof = results.(ptop).profiles;
        % Ensure we have three stored profiles
        xB = S.baseline.profile.x;
        yB = S.baseline.profile.zr;
        % Find −20% and +20%
        yM = []; yP = [];
        for k=1:numel(prof)
            if abs(prof(k).scale-0.80)<1e-9, yM = prof(k).zr; end
            if abs(prof(k).scale-1.20)<1e-9, yP = prof(k).zr; end
        end
        if ~isempty(yM) && ~isempty(yP)
            f3 = mkfig(sprintf('Overlay — most sensitive: %s', upper(ptop)), [170 120 1180 760]);
            tlo3 = tiledlayout(f3,2,1,'TileSpacing','compact','Padding','compact');

            nexttile(tlo3,1);
            plot(x_true, zr_true, 'k', 'LineWidth', 1.2); hold on;
            plot(xB, yB, 'b', 'LineWidth', 1.2);
            plot(xB, yM, 'r--', 'LineWidth', 1.0);
            plot(xB, yP, 'm--', 'LineWidth', 1.0);
            grid on; box on; formatAxes();
            title(sprintf('Truth & Recon — %s sweep (−20%%, base, +20%%)', upper(ptop)));
            xlabel('Distance (m)'); ylabel('Elevation (m)');
            legend('Truth','Base','-20%','+20%','Location','best');

            nexttile(tlo3,2);
            eB = yB - safe_interp(x_true, zr_true, xB);
            eM = yM - safe_interp(x_true, zr_true, xB);
            eP = yP - safe_interp(x_true, zr_true, xB);
            plot(xB, eB, 'b', 'LineWidth', 1.0); hold on;
            plot(xB, eM, 'r--', 'LineWidth', 1.0);
            plot(xB, eP, 'm--', 'LineWidth', 1.0);
            grid on; box on; formatAxes();
            title('Error traces  (Recon - Truth)');
            xlabel('Distance (m)'); ylabel('Error (m)'); yline(0,'k:');
            legend('Base','-20%','+20%','Location','best');

            if export_png, exportgraphics(f3, fullfile(fig_dir,'22_overlay_most_sensitive.png'), 'Resolution', 180); end
        end
    end
end

%% =============================== FUNCTIONS ===============================
function val = pickval(x, fallback)
    if isempty(x), val = fallback; else, val = x; end
end

function rec = reconstruct_once(t, azs, fs, V, x_t, x_grid, P, ...
                                hp_fc, hp_or, useSG, win_sec, poly, ...
                                lp_order, restore_mean)
    % End-to-end reconstruction for given parameter set P.
    % 1) preprocess (de-bias, integrate, high-pass)
    % 2) solve for z_u
    % 3) differentiate to z_u̇, z_ü (smoothed)
    % 4) solve for z_r(t)
    % 5) map to spatial and low-pass with n_c = f_wh/3V

    azs = azs(:) - mean(azs,'omitnan');
    t   = t(:); x_t = x_t(:); x_grid = x_grid(:);

    % Integrate
    vzs_raw = cumtrapz(t, azs);
    zss_raw = cumtrapz(t, vzs_raw);

    % High-pass drift control
    [zss, vzs] = hp_disp(zss_raw, vzs_raw, fs, hp_fc, hp_or);

    % Sprung internal force
    fs_in = P.ms*azs + P.cs*vzs + P.ks*zss;

    % Solve Cs*z_u̇ + Ks*z_u = f_s
    cs_eff = max(P.cs, 1e-6); % Add protection
    ct_eff = max(P.ct, 1e-6); % Add protection

    A1 = -(P.ks/cs_eff); B1 = (1/cs_eff); C1 = 1; D1 = 0;
    sys1 = ss(A1,B1,C1,D1);
    zu   = lsim(sys1, fs_in, t);

    % Derivatives of z_u
    vzu = smooth_derivative(zu, t, useSG, win_sec, poly);
    azu = smooth_derivative(vzu, t, useSG, win_sec, poly);

    % Solve Ct*z_ṙ + Kt*z_r = Mu*z_ü + (Cs+Ct)z_u̇ + (Ks+Kt)z_u - (Cs*z_ṡ + Ks*z_s)
    rhs = P.mu*azu + (cs_eff+ct_eff)*vzu + (P.ks+P.kt)*zu - (P.cs*vzs + P.ks*zss);
    A2 = -(P.kt/ct_eff); B2 = (1/ct_eff); C2 = 1; D2 = 0;
    sys2 = ss(A2,B2,C2,D2);
    zr_t = lsim(sys2, rhs, t);

    % Spatial mapping & low-pass
    % Compute cutoff nc from wheel-hop/3
    f_wh = (1/(2*pi))*sqrt(P.kt/P.mu);          % Hz
    nc   = f_wh / 3 / max(V,eps);               % cycles/m

    zr_x_raw = interp1(x_t, zr_t, x_grid, 'linear', 'extrap'); zr_x_raw(isnan(zr_x_raw)) = 0;

    dx         = mean(diff(x_grid),'omitnan');
    fs_spatial = 1/max(dx,eps);
    nyq_sp     = fs_spatial/2;

    Wn = clamp01(nc/nyq_sp);
    [b_lp,a_lp] = butter(max(1,lp_order), Wn, 'low');
    zr_x = filtfilt(b_lp,a_lp, zr_x_raw);

    if restore_mean
        zr_x = zr_x + (mean(zr_x_raw,'omitnan') - mean(zr_x,'omitnan'));
    end

    rec = struct();
    rec.zr_x    = zr_x(:);
    rec.zr_x_raw= zr_x_raw(:);
    rec.nc      = nc;
    rec.fs_spatial = fs_spatial;
end

function [M, bands_used] = compare_metrics(x_rec, zr_rec, x_true, zr_true, fs_spatial, remove_dc_for_psd)
    % Interpolate reconstruction onto truth grid and compute metrics.
    xmin = max(min(x_true), min(x_rec));
    xmax = min(max(x_true), max(x_rec));
    mT = (x_true>=xmin) & (x_true<=xmax);
    mR = (x_rec >=xmin) & (x_rec <=xmax);
    xT = x_true(mT); yT = zr_true(mT);
    xR = x_rec(mR);  yR = zr_rec(mR);

    y_rec_on_true = safe_interp(xR, yR, xT);
    e = y_rec_on_true - yT;

    rmse = sqrt(mean(e.^2, 'omitnan'));
    mae  = mean(abs(e), 'omitnan');
    bias = mean(e, 'omitnan');
    r    = robust_corr(yT, y_rec_on_true);

    % Spectral RSE (overall)
    [Ptrue, Ntrue] = welch_psd_spatial(yT, fs_spatial, [], 0.5, 0.5, remove_dc_for_psd);
    [Prec,  Nrec ] = welch_psd_spatial(y_rec_on_true, fs_spatial, [], 0.5, 0.5, remove_dc_for_psd);
    n_common = linspace(0, min([Ntrue(end), Nrec(end), fs_spatial/2]), max(numel(Ntrue),numel(Nrec)))';
    Ptrue_c = safe_interp(Ntrue, Ptrue, n_common);
    Prec_c  = safe_interp(Nrec,  Prec,  n_common);
    rse_all = sqrt( sum((Prec_c - Ptrue_c).^2) / (sum(Ptrue_c.^2) + eps) );

    M = struct();
    M.pointwise.bias_m = bias;
    M.pointwise.rmse_m = rmse;
    M.pointwise.mae_m  = mae;
    M.pointwise.pearson_r = r;
    M.spectral.rse_all = rse_all;

    bands_used = []; % (kept for future expansion)
end

function [Pxx, N, cfg_out] = welch_psd_spatial(x, fs_spatial, nfft_user, win_frac, ovlp_frac, remove_dc)
    x = x(:);
    if remove_dc
        x = x - mean(x,'omitnan');
    end
    Npts = numel(x);
    if isempty(nfft_user)
        nfft = max(1024, 2^nextpow2(min(8192, Npts)));
    else
        nfft = max(256, nfft_user);
    end
    wlen  = max(64, round(win_frac * nfft));
    wlen  = min(wlen, Npts);
    if mod(wlen,2)==0, wlen = wlen-1; end
    nover = max(0, round(ovlp_frac * wlen));
    win   = hamming(wlen);
    [Pxx, N] = pwelch(x, win, nover, nfft, fs_spatial, 'onesided');
    cfg_out = struct('nfft',nfft,'wlen',wlen,'nover',nover,'remove_dc',logical(remove_dc));
end

function [z_hp, vz_hp] = hp_disp(z_raw, vz_raw, fs, fc, order)
    if isempty(fc) || fc<=0
        z_hp  = z_raw;
        vz_hp = vz_raw;
        return;
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
        r = NaN; return;
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

function write_sensitivity_csv(filename, results)
    % Flatten results struct into CSV rows: param,scale,rmse,mae,bias,r,rse_all
    fid = fopen(filename,'w');
    assert(fid>0, 'Cannot open %s for writing.', filename);
    fprintf(fid,'param,scale,rmse,mae,bias,r,rse_all\n');
    pnames = fieldnames(results);
    for i=1:numel(pnames)
        p = pnames{i};
        sc = results.(p).scales;
        rm = results.(p).rmse;
        ma = results.(p).mae;
        bi = results.(p).bias;
        rr = results.(p).r;
        rs = results.(p).rse_all;
        for k=1:numel(sc)
            fprintf(fid,'%s,%.6f,%.9e,%.9e,%.9e,%.6f,%.6f\n', p, sc(k), rm(k), ma(k), bi(k), rr(k), rs(k));
        end
    end
    fclose(fid);
end

function barh_cool(labels, delta)
    % Tornado-style horizontal bar chart for ΔRMSE at -20% and +20%
    % delta: [N x 2] matrix, columns = [Δ at -20%, Δ at +20%]
    N = numel(labels);
    y = 1:N;
    dneg = delta(:,1);
    dpos = delta(:,2);

    hold on; grid on; box on;
    bh1 = barh(y, dneg, 0.4, 'FaceAlpha',0.85);   % -20%
    bh2 = barh(y, dpos, 0.4, 'FaceAlpha',0.85);
    bh1.FaceColor = [0.85 0.45 0.45];
    bh2.FaceColor = [0.45 0.65 0.95];

    set(gca,'YTick', y, 'YTickLabel', labels, 'YDir','reverse');
    xlabel('\Delta RMSE (m) relative to baseline');
    legend('-20%','+20%','Location','best');
    title('Tornado chart: sensitivity at \pm20%');
end