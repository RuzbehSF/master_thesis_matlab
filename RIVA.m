
%% RIVA Code: Time-Domain Profile Reconstruction Using PID Inversion Loop
% -------------------------------------------------------------------------
% PURPOSE:
% Implements a time-domain road profile reconstruction (Path E) using a
% Proportional-Integral-Derivative (PID) controller. This script is an
% adaptation of the RIVA_insilico workflow, redesigned to integrate with
% the 'recon_cfg.mat' data standard.
%
% CORE LOGIC (How it works):
% 1.  Loads the configuration file ('recon_cfg.mat') containing the sprung
%     acceleration signal a_s(t) [m/s^2], sampling rate fs, and speed V.
%
% 2.  Initializes a 'sys' structure for the RIVA functions.
%     - NOTE: This logic uses its own hard-coded lumped parameters (K1, K2,
%       C, U) and IGNORES the vehicle parameters (ms, mu, ks, etc.)
%       from the cfg file, as its internal math is different.
%     - The target acceleration 'sys.accm' is set from the input 'az_s',
%       converting units from m/s^2 to mm/s^2 for RIVA's consistency.
%
% 3.  Uses an optimization algorithm (fminsearch) to find optimal PID
%     gains (Kp, Ki, Kd) by calling 'PID_profile_inv_tvar'. This
%     function simulates the quarter-car model and uses the PID controller
%     to generate a profile 'Zp' (in mm) that minimizes the error between
%     its own acceleration and the target acceleration.
%
% 4.  The final optimal gains (xmin) are used to generate the final
%     reconstructed profile 'pinv' (in mm) and the corresponding
%     acceleration 'accinv' (in mm/s^2) via 'PID_profile_out_tvar'.
%
% 5.  The reconstructed profile 'pinv' is converted from mm back to meters.
%
% 6.  The result is mapped to a spatial grid 'x_grid' based on the
%     'spatial_grid_policy' (just like Path B).
%
% 7.  Saves a standard spatial struct 'sp' (x, zr_x_raw, zr_x_filt, dx,
%     fs_spatial, metadata) to the Path E output folder. The metadata
%     now includes the optimal Kp, Ki, and Kd values.
%
% 8.  (Optional) Generates publication-ready plots: Acceleration
%     comparison, PID gain sensitivity, and the standard spatial
%     comparison vs. ground truth (if available).
%
% INPUTS:
% - '00_Outputs/02_ReconstructionSetup/Config/recon_cfg.mat':
%     - 'cfg.t': Time vector [s]
%     - 'cfg.fs': Sampling rate [Hz]
%     - 'cfg.V': Vehicle speed [m/s]
%     - 'cfg.az_s' (or similar): Sprung acceleration [m/s^2]
%     - 'cfg.ground' (optional): Ground truth for plotting.
%
% OUTPUTS:
% - '03E_PID Inversion/recon_spatial.mat' : Spatial profile struct 'sp' for Path E.
% - (Optional) '03E_PID Inversion/figs/*' : Saved figures (PNG/FIG/PDF/TIFF).
% -------------------------------------------------------------------------

close all; clear; clc;
fprintf('\n=== PATH E: Time-Domain PID Inversion ===\n');

%% ---------------------------- HYPERPARAMETERS -----------------------------
% Initial PID coefficients guess
Kp0 = 0;
Ki0 = 0;
Kd0 = 0;
x0  = [Kp0 Ki0 Kd0];

% Optimization options
fmin_options = optimset('MaxIter',1e6, 'MaxFunEvals',1e9, ...
                        'TolFun',1e-4, 'TolX',1e-4, ...
                        'Display','iter', ...          % built-in text output
                        'OutputFcn',@nm_outfun);       % custom visual + verbose

% --- Workflow Parameters (from Path B) ---
spatial_grid_policy   = 'ground_truth'; % {'ground_truth','uniform'}
publish_to_shared     = true;

% ---- Plot controls ----
make_plots            = true;
save_plots            = true;      % write PNG + FIG under out_dir_pathE/figs
fig_format            = 'png';     % {'png','tiff','pdf','png'}
fig_dpi               = 300;       % export resolution
use_latex_labels      = false;     % set true if you want LaTeX (requires proper setup)

% --- Output Directories ---
out_dir_pathE = fullfile('00_Outputs', '03_ReconstructionCore', 'PathE_TimeDomainPID');
fig_dir       = fullfile(out_dir_pathE,'figs');
if ~exist(out_dir_pathE,'dir')
    mkdir(out_dir_pathE);
end
if make_plots && save_plots && ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% -------------------------- INPUT AUTO-DISCOVERY --------------------------
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};
cfg_path = '';

for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file')
        cfg_path = cfg_candidates{k};
        break;
    end
end

assert(~isempty(cfg_path), 'PATH E: recon_cfg.mat not found. Run prepare_reconstruction_inputs.m first.');
S = load(cfg_path);  assert(isfield(S,'cfg'), 'PATH E: recon_cfg.mat missing cfg struct.');
cfg = S.cfg;
t   = cfg.t(:);
fs  = cfg.fs;
V   = cfg.V;

% Find acceleration field
accel_field_order = {'az_s','a_z_s','a_s','azs','a_sprung'};
az_s = [];

for f = accel_field_order
    if isfield(cfg, f{1})
        az_s = cfg.(f{1})(:);
        break;
    end
end

assert(~isempty(az_s),'PATH E: Could not locate sprung acceleration in cfg (tried %s).', strjoin(accel_field_order,', '));

% Find native x-grid (for resampling later)
if isfield(cfg,'x_t') && ~isempty(cfg.x_t)
    x_t_input = cfg.x_t(:);
else
    x_t_input = (t - t(1)) * V;
end

N = numel(t);
fprintf('Inputs: N=%d, fs=%.3f Hz, V=%.3f m/s\n', N, fs, V);

%% ------------------- 1) SETUP 'sys' STRUCT FOR RIVA --------------------
% This bridges the 'cfg' inputs to the 'sys' struct required by RIVA's
% local functions.
fprintf('Setting up RIVA (PID) system struct...\n');
sys = struct(); % <-- Initialize sys as a new struct
sys.dx = V / fs;
sys.dt = 1 / fs;
sys.v  = V;
sys.time = t; % Use time vector from cfg

% Load vehicle parameters from cfg
sys.ms = cfg.ms;
sys.mu = cfg.mu;
sys.ks = cfg.ks;
sys.cs = cfg.cs;
sys.kt = cfg.kt;
sys.ct = cfg.ct;
fprintf('  Using vehicle parameters (ms, mu, ks, etc.) from recon_cfg.mat.\n');

% Set the target acceleration from cfg (units are m/s^2)
sys.accm = az_s; 
fprintf('  Using target acceleration in m/s^2.\n');

% Initialize response vectors
sys.Zu0  = zeros(N,1);            
sys.Zs0  = sys.Zu0; 
sys.Zp0  = sys.Zu0; 
sys.acc0 = sys.Zu0;

% Set solver type
sys.solver = 'nms';

%% -------------------- 2) RUN PID OPTIMIZATION (INVERSION) -----------------
% Run optimization algorithm
fprintf('Running fminsearch to find optimal PID gains...\n');
[xmin,fval,exitflag,output] = fminsearch(@(x)PID_profile_inv_tvar(x,sys),x0,fmin_options); % Call to local function
fprintf('Optimization complete.\n');
fprintf('  Optimal Kp: %.4f\n', xmin(1));
fprintf('  Optimal Ki: %.4f\n', xmin(2));
fprintf('  Optimal Kd: %.4f\n', xmin(3));
fprintf('  Final residual (fval): %.4g\n', fval);

%% --------------------- 3) RECONSTRUCT FINAL PROFILE ---------------------
% Get optimal profile using the gains found by the optimizer
% 'pinv' will be in [m]
% 'accinv' will be in [m/s^2]
[o,pinv,accinv] = PID_profile_out_tvar(xmin,sys); % Call to local function
fprintf('Generated final profile using optimal gains.\n');

%% --------------------- 4) MAP TO SPATIAL STRUCT & SAVE ------------------
% RIVA output 'pinv' is now in [m], no conversion needed
zr_t_meters = pinv; % [m]

% This is the "native" distance grid for the time-domain result
if numel(x_t_input) ~= N
    x_t_native = linspace(0, (N-1)*sys.dx, N)';
    warning('PATH E: Cfg x_t and t length mismatch. Using uniform dx grid.');
else
    x_t_native = x_t_input;
end

% Resample the result onto the chosen spatial grid (from Path B logic)
switch lower(spatial_grid_policy)
    case 'ground_truth'
        if isfield(cfg,'ground') && isfield(cfg.ground,'x')
            x_grid = cfg.ground.x(:);
        else
            warning('PATH E: ground truth grid missing; falling back to UNIFORM.');
            x_grid = x_t_native;
        end
    case 'uniform'
        x_grid = x_t_native;
    otherwise
        error('PATH E: Unknown spatial_grid_policy: %s', spatial_grid_policy);
end

zr_x_raw  = interp1(x_t_native, zr_t_meters, x_grid, 'linear', 'extrap');
zr_x_raw(~isfinite(zr_x_raw)) = 0;

dx         = mean(diff(x_grid),'omitnan');
fs_spatial = 1 / max(dx, eps);

% --- Spatial high-pass filter to remove wavelengths > 50 m (as in RIVA paper) ---
lambda_c      = 200;              % [m] cutoff wavelength
fc_hp         = 1 / lambda_c;    % [cycles/m]
nyq_sp        = fs_spatial / 2;  % Nyquist spatial frequency [cycles/m]
Wn_hp         = fc_hp / max(nyq_sp, eps);  % normalized cutoff [0,1]
hp_order      = NaN;
filter_applied = false;

if isfinite(Wn_hp) && Wn_hp > 0 && Wn_hp < 1
    hp_order = 4; % 4th-order Butterworth high-pass
    [b_hp, a_hp] = butter(hp_order, Wn_hp, 'high');
    zr_x_filt = filtfilt(b_hp, a_hp, zr_x_raw);
    filter_applied = true;
else
    % If dx is weird or cutoff not resolvable, fall back to raw
    zr_x_filt = zr_x_raw;
end

% --- Create the standard 'sp' struct ---
sp = struct();
sp.x              = x_grid;
sp.zr_x_raw       = zr_x_raw;
sp.zr_x_filt      = zr_x_filt;
sp.dx             = dx;
sp.fs_spatial     = fs_spatial;
sp.nc_spatial     = 1 / lambda_c;   % cutoff spatial frequency [cycles/m]
sp.lp_order       = hp_order;       % reuse field name for filter order
sp.filter_applied = filter_applied;

sp.meta = struct( ...
    'recon_method',         'Path E: Time-Domain PID Inversion', ...
    'optimal_Kp',            xmin(1), ...
    'optimal_Ki',            xmin(2), ...
    'optimal_Kd',            xmin(3), ...
    'final_residual',        fval, ...
    'grid_policy',           spatial_grid_policy, ...
    'fs',                    fs, ...
    'V',                     V ...
);

% Save Path E artifact
out_mat_pathE = fullfile(out_dir_pathE,'recon_spatial.mat');
save(out_mat_pathE, 'sp','-v7.3');
fprintf('PATH E: saved %s\n', out_mat_pathE);

% Optional publish for downstream compatibility
if publish_to_shared
    share_dir = fullfile('00_Outputs', '03_ReconstructionCore', 'SpatialMapping');
    if ~exist(share_dir,'dir')
        mkdir(share_dir);
    end
    copyfile(out_mat_pathE, fullfile(share_dir,'recon_spatial.mat'));
    fprintf('PATH E: also published to %s\n', share_dir);
end

%% ------------------------------ 5) PLOTS ----------------------------------
if make_plots
    % ---------- Global style (from Path B) ----------
    baseFont     = 11;
    axesLW       = 0.9;
    lineLW       = 1.6;
    useMinorGrid = true;

    if use_latex_labels
        set(0,'defaultTextInterpreter','latex');
        set(0,'defaultLegendInterpreter','latex');
        set(0,'defaultAxesTickLabelInterpreter','latex');
    end
    set(0,'defaultAxesFontName','Helvetica','defaultAxesFontSize',baseFont);
    set(0,'defaultLineLineWidth',lineLW);
    set(0,'DefaultAxesBox','on');

    % Color palette
    C.black = [0.10 0.10 0.10];
    C.blue  = [0.00 0.45 0.74];
    C.red   = [0.85 0.33 0.10];
    C.gray  = [0.60 0.60 0.60];
    C.green = [0.20 0.62 0.20];
    C.purp  = [0.49 0.18 0.56];

    %% ===== Fig 1: Acceleration + Residual (RIVA specific) =====
    e_acc = az_s(:) - accinv(:);
    rmse_acc = sqrt(mean(e_acc.^2,'omitnan'));
    mae_acc  = mean(abs(e_acc),'omitnan');

    fig1 = figure('Name','Path E — Acceleration & Residual','Color','w','Units','pixels');
    fig1.Position(3:4) = [980 520];
    tlo1 = tiledlayout(fig1,2,1,'TileSpacing','compact','Padding','compact');

    % Top: acceleration comparison
    ax1 = nexttile(tlo1,1); hold(ax1,'on');
    plot(ax1, t, az_s, '-',  'Color', C.black, 'DisplayName','Target a_s(t)');
    plot(ax1, t, accinv, '--','Color', C.blue,  'DisplayName','Recon. a_s(t)');
    applyTheme(ax1, useMinorGrid);
    xlabel(ax1,'Time [s]');
    ylabel(ax1,'Acceleration [m/s^2]');
    title(ax1,sprintf('Target vs reconstructed acceleration (RMSE = %.3g m/s^2)', rmse_acc));
    niceLegend(ax1,'best');

    % Bottom: acceleration residual
    ax1b = nexttile(tlo1,2); hold(ax1b,'on');
    plot(ax1b, t, e_acc, '-', 'Color', C.red, 'DisplayName','Residual a_s^{meas} - a_s^{recon}');
    yline(ax1b,0,':','Color',C.gray,'HandleVisibility','off');
    applyTheme(ax1b, useMinorGrid);
    xlabel(ax1b,'Time [s]');
    ylabel(ax1b,'Residual [m/s^2]');
    title(ax1b,sprintf('Acceleration residual (RMSE = %.3g, MAE = %.3g)', rmse_acc, mae_acc));
    niceLegend(ax1b,'best');

    linkaxes([ax1,ax1b],'x');
    maybe_export_pub(fig1, fig_dir, 'E01_Acceleration_Residual', fig_format, fig_dpi, save_plots);

    %% ===== Fig 2: PID Sensitivity (normalised objective) =====
    fig2 = figure('Name','Path E — PID Sensitivity','Color','w','Units','pixels');
    fig2.Position(3:4) = [980 400];
    tlo2 = tiledlayout(fig2,1,1,'TileSpacing','compact','Padding','compact');
    ax2 = nexttile(tlo2,1); hold(ax2,'on');

    % --- Sensitivity calculation (from original RIVA, but normalised) ---
    ulim   = 1; 
    llim   = -1; 
    pts    = 101; % Reduced points for faster plotting
    dist   = linspace(llim,ulim,pts);

    Kpe    = xmin(1) + xmin(1).*dist;
    Kie    = xmin(2) + xmin(2).*dist;
    Kde    = xmin(3) + xmin(3).*dist;

    op     = zeros(length(dist),1); 
    oi     = op; 
    od     = op;

    % Baseline objective at optimal gains (xmin)
    psi0 = PID_profile_out_tvar(xmin,sys);

    for i = 1:length(Kpe)
        x00   = [Kpe(i) xmin(2) xmin(3)];
        op(i) = PID_profile_out_tvar(x00,sys);
    end
    for i = 1:length(Kie)
        x00   = [xmin(1) Kie(i) xmin(3)];
        oi(i) = PID_profile_out_tvar(x00,sys);
    end
    for i = 1:length(Kde)
        x00   = [xmin(1) xmin(2) Kde(i)];
        od(i) = PID_profile_out_tvar(x00,sys);
    end

    % Stack, remove NaNs, normalise by psi0
    opm = [Kpe' op]; opm(any(isnan(opm),2),:) = [];
    oim = [Kie' oi]; oim(any(isnan(oim),2),:) = [];
    odm = [Kde' od]; odm(any(isnan(odm),2),:) = [];

    opm(:,2) = opm(:,2) ./ max(eps,psi0);
    oim(:,2) = oim(:,2) ./ max(eps,psi0);
    odm(:,2) = odm(:,2) ./ max(eps,psi0);

    % Clip insane values (loss of stability) so they don’t blow up the plot
    psi_min = 1e-2;
    psi_max = 1e3;

    maskP = isfinite(opm(:,2)) & opm(:,2) >= psi_min & opm(:,2) <= psi_max;
    maskI = isfinite(oim(:,2)) & oim(:,2) >= psi_min & oim(:,2) <= psi_max;
    maskD = isfinite(odm(:,2)) & odm(:,2) >= psi_min & odm(:,2) <= psi_max;

    % Plot only the valid region; disappearance of the curve ~= "unstable/outside range"
    plot(ax2, opm(maskP,1)./max(eps,xmin(1)), opm(maskP,2), '--',...
        'Color',C.blue,  'LineWidth',1.5, 'DisplayName','K_P^{FS}/K_P^0');
    plot(ax2, oim(maskI,1)./max(eps,xmin(2)), oim(maskI,2), '-.',...
        'Color',C.red,   'LineWidth',1.5, 'DisplayName','K_I^{FS}/K_I^0');
    plot(ax2, odm(maskD,1)./max(eps,xmin(3)), odm(maskD,2), ':',...
        'Color',C.black, 'LineWidth',1.5, 'DisplayName','K_D^{FS}/K_D^0');

    % Mark optimal gains
    xline(ax2,1,':','Color',C.gray,'HandleVisibility','off');

    applyTheme(ax2, useMinorGrid);
    xlabel(ax2,'K^{FS}/K^0  (fractional change from optimal gain)');
    ylabel(ax2,'Normalised objective |\psi| / |\psi_0|');
    title(ax2,'PID gain sensitivity (normalised objective)');
    niceLegend(ax2,'best');
    set(ax2,'YScale','log');
    xlim(ax2, [1+llim 1+ulim]);
    ylim(ax2, [psi_min psi_max]);

    maybe_export_pub(fig2, fig_dir, 'E02_PID_Sensitivity', fig_format, fig_dpi, save_plots);

    %% ===== Fig 3: Spatial profile vs ground + error panel =====
    if isfield(cfg,'ground') && isfield(cfg.ground,'x') && isfield(cfg.ground,'zr')
        xg = cfg.ground.x(:);
        zg = cfg.ground.zr(:);

        [x_common, zr_est_resampled, zg_resampled] = ...
            align_on_common_grid(sp.x, sp.zr_x_filt, xg, zg);

        err = zr_est_resampled - zg_resampled;

        % Metrics
        rmse = sqrt(mean(err.^2,'omitnan'));
        mae  = mean(abs(err),'omitnan');
        vp   = isfinite(zr_est_resampled) & isfinite(zg_resampled);

        if nnz(vp) >= 2 && std(zr_est_resampled(vp)) > 0 && std(zg_resampled(vp)) > 0
            rho = corr(zr_est_resampled(vp), zg_resampled(vp));
        else
            rho = NaN;
        end

        emax = robust_max(abs(err));
        emax = max(emax, eps);
        elims = 1.05*[-emax, emax];

        fig3 = figure('Name','Path E — Spatial Comparison','Color','w','Units','pixels');
        fig3.Position(3:4) = [980 620];
        tlo3 = tiledlayout(fig3,2,1,'TileSpacing','compact','Padding','compact');

        % Top: ground vs recon
        ax5 = nexttile(tlo3,1); hold(ax5,'on');
        plot(ax5, xg, zg, '-', 'Color', C.black, 'DisplayName','Ground truth');
        plot(ax5, sp.x, sp.zr_x_filt, '-', 'Color', C.blue,  'DisplayName','PID inversion (Path E)');
        applyTheme(ax5, useMinorGrid);
        ylabel(ax5,'Elevation z_r(x) [m]');
        title(ax5,'Spatial profile (ground vs RIVA)');
        niceLegend(ax5,'best');

        % Bottom: error
        ax6 = nexttile(tlo3,2); hold(ax6,'on');
        plot(ax6, x_common, err, '-', 'Color', C.red, 'DisplayName','Error (recon − truth)');
        yline(ax6, 0, ':', 'Color', C.gray, 'HandleVisibility','off');
        applyTheme(ax6, useMinorGrid);
        xlabel(ax6,'Distance x [m]');
        ylabel(ax6,'Error [m]');
        ylim(ax6, elims);
        yline(ax6,0,':','Color',C.gray,'HandleVisibility','off');
        title(ax6, sprintf('Reconstruction error (RMSE = %.3g m, MAE = %.3g m, \\rho = %.3f)', ...
        rmse, mae, rho));
        linkaxes([ax5,ax6],'x');
        niceLegend(ax6,'best');

        maybe_export_pub(fig3, fig_dir, 'E03_Spatial_Comparison', fig_format, fig_dpi, save_plots);
    end

    %% ===== Fig 4: Time-domain sanity checks =====
    fig4 = figure('Name','Path E — Time Signals','Color','w','Units','pixels');
    fig4.Position(3:4) = [980 620];
    tlo4 = tiledlayout(fig4,2,1,'TileSpacing','compact','Padding','compact');

    % Top: measured vs reconstructed acceleration (zoomed visual)
    ax7 = nexttile(tlo4,1); hold(ax7,'on');
    plot(ax7, t, az_s, '-',  'Color', C.black, 'DisplayName','a_s^{meas}(t)');
    plot(ax7, t, accinv, '--','Color', C.blue,  'DisplayName','a_s^{recon}(t)');
    applyTheme(ax7, useMinorGrid);
    ylabel(ax7,'Acceleration [m/s^2]');
    title(ax7,'Sprung acceleration (time domain, overlay)');
    niceLegend(ax7,'best');

    % Bottom: reconstructed road in spatial domain (native grid)
    ax8 = nexttile(tlo4,2); hold(ax8,'on');
    plot(ax8, x_t_native, zr_t_meters, 'Color', C.green, 'DisplayName','z_r(t) mapped to x');
    applyTheme(ax8, useMinorGrid);
    xlabel(ax8,'Distance x = V(t - t_0) [m]');
    ylabel(ax8,'Reconstructed z_r [m]');
    title(ax8,'Reconstructed road before spatial resampling');
    niceLegend(ax8,'best');

    maybe_export_pub(fig4, fig_dir, 'E04_TimeSignals', fig_format, fig_dpi, save_plots);

end
fprintf('=== PATH E: done ===\n\n');

%% --------------------- HELPERS (LOCAL FUNCTIONS) ------------------------

function stop = nm_outfun(x, optimValues, state)
% Output function for fminsearch (Nelder–Mead) to provide:
%   - Command-window verbose logging
%   - Live plots of:
%       * Objective value vs iteration
%       * Kp vs iteration
%       * Ki vs iteration
%       * Kd vs iteration

    persistent figNm ax1 ax2 ax3 ax4 ...
               lineF lineKp lineKi lineKd ...
               itHist fHist kpHist kiHist kdHist

    stop = false;  % never stop prematurely

    switch state
        case 'init'
            % Initialise history
            itHist = [];
            fHist  = [];
            kpHist = [];
            kiHist = [];
            kdHist = [];

            % Create figure for visualisation
            figNm = figure('Name','Nelder-Mead optimisation progress', ...
                           'NumberTitle','off', ...
                           'Color','w');

            % Layout: 4x1 subplots
            % 1) Objective
            ax1 = subplot(4,1,1,'Parent',figNm);
            hold(ax1,'on');
            lineF = plot(ax1, nan, nan, '-');
            xlabel(ax1,'Iteration');
            ylabel(ax1,'Objective f');
            title(ax1,'Objective value vs iteration');
            grid(ax1,'on');

            % 2) Kp
            ax2 = subplot(4,1,2,'Parent',figNm);
            hold(ax2,'on');
            lineKp = plot(ax2, nan, nan, '-');
            xlabel(ax2,'Iteration');
            ylabel(ax2,'K_p');
            title(ax2,'K_p vs iteration');
            grid(ax2,'on');

            % 3) Ki
            ax3 = subplot(4,1,3,'Parent',figNm);
            hold(ax3,'on');
            lineKi = plot(ax3, nan, nan, '-');
            xlabel(ax3,'Iteration');
            ylabel(ax3,'K_i');
            title(ax3,'K_i vs iteration');
            grid(ax3,'on');

            % 4) Kd
            ax4 = subplot(4,1,4,'Parent',figNm);
            hold(ax4,'on');
            lineKd = plot(ax4, nan, nan, '-');
            xlabel(ax4,'Iteration');
            ylabel(ax4,'K_d');
            title(ax4,'K_d vs iteration');
            grid(ax4,'on');

        case 'iter'
            k  = optimValues.iteration;
            f  = optimValues.fval;

            % Append to history
            itHist(end+1) = k;
            fHist(end+1)  = f;
            kpHist(end+1) = x(1);
            kiHist(end+1) = x(2);
            kdHist(end+1) = x(3);

            % Update plots if figure is still valid
            if ~isempty(figNm) && isvalid(figNm)
                set(lineF,  'XData', itHist, 'YData', fHist);
                set(lineKp, 'XData', itHist, 'YData', kpHist);
                set(lineKi, 'XData', itHist, 'YData', kiHist);
                set(lineKd, 'XData', itHist, 'YData', kdHist);

                drawnow limitrate;
            end

            % Command window verbose
            fprintf('  [Nelder-Mead] iter %4d | f = %.4g | Kp = %.4g | Ki = %.4g | Kd = %.4g\n', ...
                    k, f, x(1), x(2), x(3));

        case 'done'
            % Final summary
            fprintf('  [Nelder-Mead] finished at iter %4d | f = %.4g | Kp = %.4g | Ki = %.4g | Kd = %.4g\n', ...
                    optimValues.iteration, optimValues.fval, x(1), x(2), x(3));

            % Keep figure open, just clear persistent vars
            clear figNm ax1 ax2 ax3 ax4 ...
                  lineF lineKp lineKi lineKd ...
                  itHist fHist kpHist kiHist kdHist;
    end
end

function o = PID_profile_inv_tvar(X,sys)
% Inversion function for fminsearch: Calculates profile & total error
% --- CORRECTED LOGIC for discrete-time state-space ---

% Identify variable parameters
Kp = X(1);
Ki = X(2);
Kd = X(3);

% Load parameters from sys struct
ms = sys.ms; mu = sys.mu;
ks = sys.ks; cs = sys.cs;
kt = sys.kt; ct = sys.ct;

% Get time vector and target acceleration
time = sys.time;
Zsm_dotdot_target = sys.accm; % Target acceleration [m/s^2]
n = length(Zsm_dotdot_target);
dt = time(2) - time(1); % Assume constant dt

% --- Define the State-Space Model ---
% x = [zs; zs_dot; zu; zu_dot]
% u = [Zp] (We assume Zp_dot = 0)
A = [0 1 0 0;
     -ks/ms -cs/ms  ks/ms   cs/ms;
      0 0 0 1;
      ks/mu cs/mu -(ks+kt)/mu -(cs+ct)/mu];
 
B_full = [0 0; 0 0; 0 0; kt/mu ct/mu]; % B for u = [Zp; Zp_dot]
B = B_full(:, 1); % B for u = [Zp] only

% Output matrix (C) for y = zs_dotdot
C = [-ks/ms -cs/ms ks/ms cs/ms];
D = 0; % D for u = [Zp]

% Create a discrete-time state-space system
sys_d = c2d(ss(A,B,C,D), dt, 'zoh');
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

% --- Run the PID simulation ---
Zp_pid = zeros(n, 1); % The road profile we are solving for
Zs_dotdot_calc = zeros(n, 1); % The calculated acceleration
ep = zeros(n, 1);
ei = zeros(n, 1);
ed = zeros(n, 1);
x_state = zeros(size(Ad,1), 1); % Initial state x(1) = [0;0;0;0]

% Loop from the first time step
for i = 2:n
    % 1. Evolve the state
    % x(i) = A*x(i-1) + B*u(i-1)
    % Note: x_state holds x(i-1) at the start of the loop
    x_state_new = Ad * x_state + Bd * Zp_pid(i-1);
    
    % 2. Predict the output *without* the current input
    % This is the system's response to past events
    y_predict = Cd * x_state_new;
    
    % 3. Calculate error based on this prediction
    % e_P(i) = r(i) - y_predict(i)
    ep(i) = Zsm_dotdot_target(i) - y_predict;
    
    % 4. PID logic
    % 4a) Integral term: trapezoidal integration of e_P
    ei(i) = ei(i-1) + 0.5 * dt * (ep(i) + ep(i-1));
    
    % 4b) Derivative term: difference of measured vs model accel derivatives
    %     (uses central diff for measured accel, backward diff for model)
    if i < n
        d_m = (Zsm_dotdot_target(i+1) - Zsm_dotdot_target(i-1)) / (2*dt);
    else
        d_m = (Zsm_dotdot_target(i) - Zsm_dotdot_target(i-1)) / dt;
    end
    d_c = (y_predict - Zs_dotdot_calc(i-1)) / dt;
    ed(i) = d_m - d_c;
    
    % 5. Calculate the *current* input u(i) based on the error
    % This is the PID's corrective action
    Zp_pid(i) = Kp*ep(i) + Ki*ei(i) + Kd*ed(i);
    
    % 6. Calculate the *actual* final output (for the residual)
    % y(i) = C*x(i) + D*u(i)
    Zs_dotdot_calc(i) = Cd * x_state_new + Dd * Zp_pid(i);
    
    % 7. Store the new state for the next iteration
    x_state = x_state_new;
end

% --- Calculate final residual (L1 norm, as in RIVA paper) ---
residual = Zsm_dotdot_target - Zs_dotdot_calc;
o = sum(abs(residual));  % Sum of absolute differences
end

function [o,Zp_pid,Zs_dotdot_calc] = PID_profile_out_tvar(X,sys)
% Output function: Runs the model with optimal gains to get final profile
% Identify variable parameters
Kp = X(1);
Ki = X(2);
Kd = X(3);
% Load parameters from sys struct
ms = sys.ms; mu = sys.mu;
ks = sys.ks; cs = sys.cs;
kt = sys.kt; ct = sys.ct;
% Get time vector and target acceleration
time = sys.time;
Zsm_dotdot_target = sys.accm; % Target acceleration [m/s^2]
n = length(Zsm_dotdot_target);
dt = time(2) - time(1); % Assume constant dt
% --- Define the State-Space Model ---
A = [0 1 0 0;
     -ks/ms -cs/ms  ks/ms   cs/ms;
      0 0 0 1;
      ks/mu cs/mu -(ks+kt)/mu -(cs+ct)/mu];
B_full = [0 0; 0 0; 0 0; kt/mu ct/mu];
B = B_full(:, 1); 
C = [-ks/ms -cs/ms ks/ms cs/ms];
D = 0;
sys_d = c2d(ss(A,B,C,D), dt, 'zoh');
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;
% --- Run the PID simulation ---
Zp_pid = zeros(n, 1); 
Zs_dotdot_calc = zeros(n, 1); 
ep = zeros(n, 1);
ei = zeros(n, 1);
ed = zeros(n, 1);
x_state = zeros(size(Ad,1), 1); % Initial state x(1) = [0;0;0;0]
% Loop from the first time step
for i = 2:n
    % 1. Evolve the state
    % x(i) = A*x(i-1) + B*u(i-1)
    x_state_new = Ad * x_state + Bd * Zp_pid(i-1);
    
    % 2. Predict the output *without* the current input
    y_predict = Cd * x_state_new;
    
    % 3. Calculate error based on this prediction
    ep(i) = Zsm_dotdot_target(i) - y_predict;
    
    % 4. PID logic
    % 4a) Integral term: trapezoidal integration of e_P
    ei(i) = ei(i-1) + 0.5 * dt * (ep(i) + ep(i-1));
    
    % 4b) Derivative term: difference of measured vs model accel derivatives
    if i < n
        d_m = (Zsm_dotdot_target(i+1) - Zsm_dotdot_target(i-1)) / (2*dt);
    else
        d_m = (Zsm_dotdot_target(i) - Zsm_dotdot_target(i-1)) / dt;
    end
    d_c = (y_predict - Zs_dotdot_calc(i-1)) / dt;
    ed(i) = d_m - d_c;
    
    % 5. Calculate the *current* input u(i) based on the error
    Zp_pid(i) = Kp*ep(i) + Ki*ei(i) + Kd*ed(i);
    
    % 6. Calculate the *actual* final output (for the residual)
    Zs_dotdot_calc(i) = Cd * x_state_new + Dd * Zp_pid(i);
    
    % 7. Store the new state for the next iteration
    x_state = x_state_new;
end

% --- Calculate final residual (L1 norm, as in RIVA paper) ---
residual = Zsm_dotdot_target - Zs_dotdot_calc;
o = sum(abs(residual));
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

function maybe_export_pub(fig, outdir, basename, fmt, dpi, doSave)
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