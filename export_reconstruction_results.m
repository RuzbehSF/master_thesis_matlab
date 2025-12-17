
%% 6) Consolidate and Export All Reconstruction Results
% -------------------------------------------------------------------------
% PURPOSE:
% This script acts as a final packaging step. It gathers all data from the
% key stages of the pipeline and consolidates them into a single, comprehensive
% .mat file and several easy-to-use .csv files for external analysis.
%
% CORE LOGIC (How it works):
% 1.  Loads the three main data files from the pipeline: the configuration
%     ('recon_cfg.mat'), the time-domain results ('recon_time_domain.mat'),
%     and the final spatial-domain profile ('recon_spatial.mat').
%
% 2.  Organizes all the loaded data into a single, neatly structured package
%     ('pkg') with clear sub-fields for configuration, time-domain signals,
%     spatial-domain signals, and ground truth data.
%
% 3.  Saves this consolidated 'pkg' structure to a final .mat file, which
%     serves as a complete archive of one full reconstruction run.
%
% 4.  Creates and exports several tables as .csv files for easy access in
%     other software (e.g., Excel, Python). This includes separate CSVs for
%     all time-domain signals, the spatial-domain profiles, and the original
%     ground truth profile.
%
% 5.  Generates a 'manifest.txt' file summarizing the key settings of the
%     run, such as vehicle speed and filter parameters.
%
% INPUTS:
% - '02_Reconstruction Inputs/recon_cfg.mat'
% - '03_Time Domain Inversion/recon_time_domain.mat'
% - '04_Spatial Mapping/recon_spatial.mat'
%
% OUTPUTS:
% - '06_Exported Results/reconstructed_profile.mat': The final consolidated data package.
% - '06_Exported Results/*.csv': CSV files for time, spatial, and ground truth data.
% - '06_Exported Results/recon_manifest.txt': A summary text file.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ========================== HYPERPARAMETERS ==========================
% (1) Input discovery (new locations first, then legacy)
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};
td_candidates  = { ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'PathA_TimeDomain', 'recon_time_domain.mat') ...
};
sp_candidates  = { ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'SpatialMapping', 'recon_spatial.mat'), ...
};

% (2) This script’s output folder (per-script organization)
save_folder = fullfile('00_Outputs', '06_FinalExports', 'ConsolidatedPackage');
fig_dir     = fullfile(save_folder,'figs');

% (3) CSV export controls
csv_decimal_precision  = 9;      % number of digits after decimal
write_ground_truth_csv = true;   % also export ground-truth spatial profile if present

% (4) QA plots & export
make_plots  = true;
export_png  = true;
% =====================================================================

%% ------- Locate required inputs -------
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file')==2, cfg_path = cfg_candidates{k}; break; end
end
assert(~isempty(cfg_path), 'Missing recon_cfg.mat in expected locations.');

td_path = '';
for k = 1:numel(td_candidates)
    if exist(td_candidates{k},'file')==2, td_path = td_candidates{k}; break; end
end
have_td = ~isempty(td_path);  % tolerate missing time-domain MAT

sp_path = '';
for k = 1:numel(sp_candidates)
    if exist(sp_candidates{k},'file')==2, sp_path = sp_candidates{k}; break; end
end
assert(~isempty(sp_path), 'Missing recon_spatial.mat in expected locations.');

% Make folders
if ~exist(save_folder,'dir'), mkdir(save_folder); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% ------- Load -------
C     = load(cfg_path);
cfg   = C.cfg;
SP    = load(sp_path);
sp    = SP.sp;

if have_td
    TD  = load(td_path);    recon = TD.recon;
else
    warning('recon_time_domain.mat not found. Proceeding with spatial-only fallback (no re-run).');
    recon = []; % sentinel
end

%% ------- Time-domain signals (full if available, fallback otherwise) -------
if have_td
    % From time-domain MAT
    t   = recon.t(:);
    fs  = cfg.fs;
    ms  = cfg.ms;
    mu  = cfg.mu;
    ks  = cfg.ks;
    cs  = cfg.cs;
    kt  = cfg.kt;
    ct  = cfg.ct;

    % Signals
    azs   = recon.azs(:);
    vzs   = recon.vzs(:);
    zss   = recon.zss(:);

    zu    = recon.zu(:);
    vzu   = recon.vzu(:);
    azu   = recon.azu(:);

    zr_t  = recon.zr_t(:);
    zr_dt = recon.zr_dot(:);

    % Sprung balance f_s (completeness)
    fs_in = ms*azs + cs*vzs + ks*zss;

    % Tyre/unsprung balance
    rhs_tyre      = mu*azu + (cs+ct)*vzu + (ks+kt)*zu - (cs*vzs + ks*zss);
    lhs_tyre      = ct*zr_dt + kt*zr_t;
    resid_tyre_N  = lhs_tyre - rhs_tyre;

else
    % Build minimal time arrays from spatial data so we can still export artifacts
    [t, zr_t, zr_dt] = time_from_spatial(cfg, sp);   % helper below
    fs  = cfg.fs;

    % Placeholders for unavailable QA signals
    azs = nan(size(t));
    vzs = nan(size(t));
    zss = nan(size(t));
    zu  = nan(size(t));
    vzu = nan(size(t));
    azu = nan(size(t));
    fs_in = nan(size(t));
    lhs_tyre = nan(size(t));
    rhs_tyre = nan(size(t));
    resid_tyre_N = nan(size(t));
end

%% ------- Build consolidated package (MAT) -------
pkg = struct();
pkg.cfg = cfg;

% Time-domain package (filled with NaNs if TD missing)
pkg.time.t_s          = t;
pkg.time.zr_t_m       = zr_t;
pkg.time.zr_dot_mps   = zr_dt;
pkg.time.zs_m         = zss;
pkg.time.zs_dot_mps   = vzs;
pkg.time.azs_mps2     = azs;
pkg.time.zu_m         = zu;
pkg.time.zu_dot_mps   = vzu;
pkg.time.azu_mps2     = azu;
pkg.time.fs_in_N      = fs_in;
pkg.time.tyre_lhs_N   = lhs_tyre;
pkg.time.tyre_rhs_N   = rhs_tyre;
pkg.time.tyre_resid_N = resid_tyre_N;

% Spatial package
pkg.spatial.x_m           = sp.x(:);
pkg.spatial.zr_x_raw_m    = sp.zr_x_raw(:);
pkg.spatial.zr_x_filt_m   = sp.zr_x_filt(:);
pkg.spatial.dx_m          = sp.dx;
pkg.spatial.fs_spatial_hz = sp.fs_spatial;   % samples per meter
pkg.spatial.nc_cyc_per_m  = sp.nc_spatial;
pkg.spatial.lp_order      = sp.lp_order;

% Ground-truth (if available)
if isfield(cfg,'ground') && ~isempty(cfg.ground) && isfield(cfg.ground,'x') && ~isempty(cfg.ground.x)
    pkg.truth.x_m = cfg.ground.x(:);
    pkg.truth.zr_true_m = cfg.ground.zr(:);
else
    pkg.truth = struct();
end

save(fullfile(save_folder,'reconstructed_profile.mat'), 'pkg', '-v7.3');
fprintf('Saved: %s\n', fullfile(save_folder,'reconstructed_profile.mat'));

%% ------- Write CSVs -------
% Time-domain CSV (full if TD present, lean otherwise)
if have_td
    T_time = table(pkg.time.t_s, pkg.time.zr_t_m, pkg.time.zr_dot_mps, ...
                   pkg.time.zs_m, pkg.time.zs_dot_mps, pkg.time.azs_mps2, ...
                   pkg.time.zu_m, pkg.time.zu_dot_mps, pkg.time.azu_mps2, ...
                   pkg.time.fs_in_N, pkg.time.tyre_lhs_N, pkg.time.tyre_rhs_N, pkg.time.tyre_resid_N, ...
        'VariableNames', {'t_s','zr_t_m','zr_dot_mps','zs_m','zs_dot_mps','azs_mps2', ...
                          'zu_m','zu_dot_mps','azu_mps2','fs_in_N','tyre_lhs_N','tyre_rhs_N','tyre_resid_N'});
else
    % Lean CSV: only what we can guarantee without upstream TD file.
    T_time = table(pkg.time.t_s, pkg.time.zr_t_m, pkg.time.zr_dot_mps, ...
        'VariableNames', {'t_s','zr_t_m','zr_dot_mps'});
end

time_csv = fullfile(save_folder,'reconstructed_profile_time.csv');
writetable_precise(T_time, time_csv, csv_decimal_precision);
fprintf('Saved: %s\n', time_csv);

% Spatial-domain CSV
T_sp = table(pkg.spatial.x_m, pkg.spatial.zr_x_raw_m, pkg.spatial.zr_x_filt_m, ...
    'VariableNames', {'x_m','zr_x_raw_m','zr_x_filt_m'});

sp_csv = fullfile(save_folder,'reconstructed_profile_spatial.csv');
writetable_precise(T_sp, sp_csv, csv_decimal_precision);
fprintf('Saved: %s\n', sp_csv);

% Ground-truth CSV (optional)
if write_ground_truth_csv && isfield(pkg,'truth') && ~isempty(pkg.truth)
    if isfield(pkg.truth,'x_m') && ~isempty(pkg.truth.x_m)
        T_true = table(pkg.truth.x_m, pkg.truth.zr_true_m, ...
            'VariableNames', {'x_m','zr_true_m'});
        true_csv = fullfile(save_folder,'true_profile_spatial.csv');
        writetable_precise(T_true, true_csv, csv_decimal_precision);
        fprintf('Saved: %s\n', true_csv);
    end
end

% Manifest
manifest_path = fullfile(save_folder,'recon_manifest.txt');
fid = fopen(manifest_path,'w');
if fid>0
    % Gentle guards for optional fields
    if isfield(cfg,'param_source')
        param_source = cfg.param_source;
    else
        param_source = '(not set)';
    end
    if isfield(cfg,'accel_used')
        accel_used = cfg.accel_used;
    else
        accel_used = '(not set)';
    end
    if isfield(cfg,'V_kmh')
        V_kmh = cfg.V_kmh;
    else
        V_kmh = cfg.V * 3.6;
    end

    fprintf(fid,'Reconstruction Artifacts Manifest\n');
    fprintf(fid,'=================================\n');
    fprintf(fid,'Params source       : %s\n', param_source);
    fprintf(fid,'Accel used          : %s\n', accel_used);
    fprintf(fid,'Speed V (m/s)       : %.4f  (%.1f km/h)\n', cfg.V, V_kmh);
    fprintf(fid,'Temporal fs (Hz)    : %.2f\n', cfg.fs);
    fprintf(fid,'Spatial fs (1/m)    : %.4f\n', pkg.spatial.fs_spatial_hz);
    fprintf(fid,'LP cutoff (cyc/m)   : %.5f  (order %d)\n', pkg.spatial.nc_cyc_per_m, pkg.spatial.lp_order);
    fprintf(fid,'Saved files:\n');
    fprintf(fid,'  - reconstructed_profile.mat\n');
    fprintf(fid,'  - reconstructed_profile_time.csv\n');
    fprintf(fid,'  - reconstructed_profile_spatial.csv\n');
    if exist('true_csv','var')
        fprintf(fid,'  - true_profile_spatial.csv\n');
    end
    if ~have_td
        fprintf(fid,'NOTE: recon_time_domain.mat was missing; exported lean time CSV only.\n');
    end
    fclose(fid);
end
fprintf('Saved: %s\n', manifest_path);

%% ------- QA Figures -------
if make_plots
    % Fig 1 — Time-domain summary (full if TD present, lean otherwise)
    f1 = mkfig('Artifacts — Time-domain summary', [100 80 1100 780]);

    if have_td
        tlo1 = tiledlayout(f1,4,1,'TileSpacing','compact','Padding','compact');

        nexttile(tlo1,1);
        plot(t, pkg.time.azs_mps2, 'LineWidth', 1.1); grid on; box on; formatAxes();
        title('Sprung acceleration  a_z^s(t)'); xlabel('Time (s)'); ylabel('m/s^2');

        nexttile(tlo1,2);
        plot(t, pkg.time.zs_dot_mps, 'b', 'LineWidth', 1.0); hold on;
        plot(t, pkg.time.zs_m,      'g', 'LineWidth', 1.0);
        grid on; box on; formatAxes();
        title('Body kinematics  z_ṡ(t), z_s(t)'); xlabel('Time (s)'); ylabel('SI units');
        legend('z_ṡ (m/s)','z_s (m)','Location','best');

        nexttile(tlo1,3);
        plot(t, pkg.time.zu_m, 'm', 'LineWidth', 1.1); hold on;
        plot(t, pkg.time.zu_dot_mps, 'r', 'LineWidth', 0.9);
        plot(t, pkg.time.azu_mps2,   'k', 'LineWidth', 0.9);
        grid on; box on; formatAxes();
        title('Unsprung response  z_u(t), z_u̇(t), z_ü(t)'); xlabel('Time (s)'); ylabel('SI units');
        legend('z_u (m)','z_u̇ (m/s)','z_ü (m/s^2)','Location','best');

        nexttile(tlo1,4);
        plot(t, pkg.time.zr_t_m, 'LineWidth', 1.1); grid on; box on; formatAxes();
        title('Reconstructed road (time domain)  z_r(t)'); xlabel('Time (s)'); ylabel('m');
    else
        plot(t, pkg.time.zr_t_m, 'LineWidth', 1.1); grid on; box on; formatAxes();
        title('Reconstructed road (time domain)  z_r(t) — spatial fallback'); xlabel('Time (s)'); ylabel('m');
    end

    if export_png, exportgraphics(f1, fullfile(fig_dir,'10_artifact_time_signals.png'), 'Resolution', 180); end

    % Fig 2 — Tyre/unsprung equation balance & residual (only if TD present)
    if have_td
        f2 = mkfig('Artifacts — Tyre/unsprung balance', [130 110 1100 760]);
        tlo2 = tiledlayout(f2,2,1,'TileSpacing','compact','Padding','compact');

        nexttile(tlo2,1);
        plot(t, pkg.time.tyre_lhs_N, 'b-', 'LineWidth', 1.1); hold on;
        plot(t, pkg.time.tyre_rhs_N, 'r--', 'LineWidth', 1.1);
        grid on; box on; formatAxes();
        title('C_t z_ṙ + K_t z_r  vs  RHS'); xlabel('Time (s)'); ylabel('N');
        legend('LHS','RHS','Location','best');

        nexttile(tlo2,2);
        plot(t, pkg.time.tyre_resid_N, 'k-', 'LineWidth', 1.0); grid on; box on; formatAxes();
        title('Tyre equation residual'); xlabel('Time (s)'); ylabel('N');

        if export_png, exportgraphics(f2, fullfile(fig_dir,'11_tyre_balance_residual.png'), 'Resolution', 180); end
    end

    % Fig 3 — Spatial raw vs filtered (full & zoom)
    f3 = mkfig('Artifacts — Spatial profiles', [160 130 1100 780]);
    tlo3 = tiledlayout(f3,2,1,'TileSpacing','compact','Padding','compact');

    nexttile(tlo3,1);
    plot(pkg.spatial.x_m, pkg.spatial.zr_x_raw_m, 'Color',[0.7 0.8 1.0], 'LineWidth', 0.9); hold on;
    plot(pkg.spatial.x_m, pkg.spatial.zr_x_filt_m, 'b', 'LineWidth', 1.2);
    grid on; box on; formatAxes();
    title(sprintf('z_r(x): raw vs filtered (n_c=%.3f cyc/m, order %d)', ...
        pkg.spatial.nc_cyc_per_m, pkg.spatial.lp_order));
    xlabel('Distance (m)'); ylabel('Elevation (m)');
    legend('raw','filtered','Location','best');

    nexttile(tlo3,2);
    if numel(pkg.spatial.x_m) > 50
        [xZ, rawZ, filZ] = middle_window(pkg.spatial.x_m, pkg.spatial.zr_x_raw_m, pkg.spatial.zr_x_filt_m, 0.12);
        plot(xZ, rawZ, 'Color',[0.85 0.9 1.0], 'LineWidth', 1.0); hold on;
        plot(xZ, filZ, 'b', 'LineWidth', 1.2);
        grid on; box on; formatAxes();
        title('z_r(x): zoomed view (middle 12%)'); xlabel('Distance (m)'); ylabel('Elevation (m)');
        legend('raw','filtered','Location','best');
    else
        axis off; text(0.1,0.5,'Zoom skipped (short profile)','Units','normalized','FontSize',11);
    end

    if export_png, exportgraphics(f3, fullfile(fig_dir,'12_spatial_profiles.png'), 'Resolution', 180); end
end

%% =============================== FUNCTIONS ===============================
function writetable_precise(T, filename, digits)
    % Write table T to CSV with fixed numeric precision.
    % Ensures consistent formatting across MATLAB versions.
    if nargin < 3 || isempty(digits)
        digits = 9;
    end
    fmt = ['%.' num2str(digits) 'f'];

    % Build format string dynamically: numeric -> fmt, text -> %s
    varnames = T.Properties.VariableNames;
    ncol = numel(varnames);

    % Open file
    fid = fopen(filename, 'w');
    assert(fid>0, 'Cannot open file for writing: %s', filename);

    % Header
    for i = 1 : ncol
        fprintf(fid, '%s', varnames{i});
        if i < ncol
            fprintf(fid, ',');
        else
            fprintf(fid, '\n');
        end
    end

    % Rows
    for r=1:height(T)
        for c=1:ncol
            val = T{r,c};
            if iscell(val)
                sval = string(val);
                sval = replace(sval, ",", ";");
                fprintf(fid, '%s', sval);
            elseif isstring(val)
                sval = replace(val, ",", ";");
                fprintf(fid, '%s', sval);
            elseif isnumeric(val)
                if isscalar(val)
                    fprintf(fid, fmt, val);
                else
                    fprintf(fid, '"[');
                    v = val(:)'; % row
                    for k = 1 : numel(v)
                        fprintf(fid, fmt, v(k));
                        if k < numel(v)
                            fprintf(fid, ';');
                        end
                    end
                    fprintf(fid, ']"');
                end
            else
                fprintf(fid, '%s', char(val));
            end
            if c < ncol
                fprintf(fid, ',');
            else
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
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

function [xZ, y1Z, y2Z] = middle_window(x, y1, y2, frac)
    if nargin < 4 || isempty(frac)
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

function [t_est, zr_t_est, zr_dot_est] = time_from_spatial(cfg, sp)
% Map spatial profile to a time axis using constant speed V (no re-run).
    x   = sp.x(:);
    V   = max(cfg.V, eps);
    % Use cfg.t(1) if available, else start at 0
    if isfield(cfg,'t') && ~isempty(cfg.t)
        t0 = cfg.t(1);
    else
        t0 = 0;
    end
    t_est = (x - x(1))/V + t0;

    % Prefer filtered if present; else raw
    if isfield(sp,'zr_x_filt') && ~isempty(sp.zr_x_filt)
        zr_x = sp.zr_x_filt(:);
    else
        zr_x = sp.zr_x_raw(:);
    end
    zr_t_est  = zr_x;
    zr_dot_est = gradient_safe(zr_t_est, t_est);
end

function dy = gradient_safe(y, t)
    y = y(:); t = t(:);
    if numel(y) < 2
        dy = zeros(size(y));
        return;
    end
    dy = zeros(size(y));
    dt = diff(t);
    dt(dt==0) = eps;
    dy(1) = (y(2)-y(1))/dt(1);
    dy(2 : end-1) = (y(3:end)-y(1:end-2))./(t(3:end)-t(1:end-2));
    dy(end) = (y(end)-y(end-1))/dt(end);
end
