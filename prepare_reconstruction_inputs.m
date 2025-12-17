
%% 2) Prepare Inputs for Road Profile Reconstruction
% -------------------------------------------------------------------------
% PURPOSE:
% This script serves as the bridge between the raw vehicle simulation and
% the profile reconstruction algorithm. It loads simulation data, selects
% the core input signals, defines the parameters for the inverse model,
% and saves a consolidated configuration file for subsequent steps.
%
% CORE LOGIC (How it works):
% 1.  Loads the 'simulation_results.mat' file containing the full output
%     from the vehicle dynamics simulation.
%
% 2.  Extracts essential data, including the time vector, vehicle speed,
%     ground truth road profile, and crucially, the sprung mass (body)
%     vertical acceleration. It allows the user to select either the "clean"
%     or "measured" (noisy) acceleration as the primary input for the
%     inversion.
%
% 3.  Selects the quarter-car parameters (masses, stiffness, damping) that
%     will be assumed by the reconstruction algorithm. This is a key step,
%     as it allows the simulation to test the algorithm's robustness when
%     its assumed parameters do not perfectly match the true vehicle's.
%
% 4.  Using these chosen "inversion parameters," it calculates the two key
%     natural frequencies of the quarter-car model: the low-frequency body
%     motion (sprung mass mode) and the high-frequency wheel hop (unsprung
%     mass mode).
%
% 5.  Based on the calculated wheel-hop frequency, it determines the
%     recommended cutoff frequency for the spatial low-pass filter that will
%     be applied at the end of the reconstruction. This follows the
%     methodology of the source paper (Ngwangwa, 2020) to ensure the
%     reconstruction is limited to a physically meaningful bandwidth.
%
% 6.  All selected signals, inversion parameters, calculated frequencies, and
%     ground truth data are packaged into a single configuration struct ('cfg')
%     and saved to a .mat file.
%
% INPUTS:
% - '01_Quarter Car Modeling/simulation_results.mat': Raw data from the
%   vehicle simulation.
%
% OUTPUTS:
% - '02_Reconstruction Inputs/recon_cfg.mat': A consolidated configuration
%   file used by all subsequent reconstruction and analysis scripts.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ========================== HYPERPARAMETERS ==========================
% (1) Paths
sim_mat_path = fullfile('00_Outputs', '01_Simulation', 'SimulationData', 'simulation_results.mat');

% (2) Which acceleration to use from your results?  'meas' (preferred) or 'clean'
accel_channel = 'meas';   % 'meas' | 'clean'

% (3) Parameters to use for inversion
% 'golden_car'
% 'from_results'
% 'paper_experimental'
% 'paper_numerical'
param_source = 'golden_Car';

% (4) Displacement drift control (preview only here; actual filtering happens later)
hp_disp_fc_hz = 0.15;     % Hz
hp_disp_order = 2;

% (5) Reconstruction low-pass policy
lp_rec_rule = 'one_third_wheelhop';
lp_rec_order = 4;

% (6) Plots
fig_export_png = true;

%% ============================= OUTPUT FOLDERS =============================
save_folder = fullfile('00_Outputs', '02_ReconstructionSetup', 'Config');
if ~exist(save_folder,'dir')
    mkdir(save_folder);
end

fig_dir = fullfile(save_folder,'figs');

if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

cfg_path  = fullfile(save_folder,'recon_cfg.mat');
diag_txt  = fullfile(save_folder,'recon_init_diagnostics.txt');

%% ============================== LOAD INPUTS ===============================
assert(exist(sim_mat_path,'file')==2, 'File not found: %s', sim_mat_path);
L = load(sim_mat_path);
assert(isfield(L,'results'), 'Expected struct "results" in %s', sim_mat_path);
R = L.results;

% Time and sampling
t  = R.time(:);
fs = R.meta.fs;
dt = 1/fs;

% Speed & distance
V_kmh = R.meta.V_kmh;
V     = V_kmh/3.6;
x_t   = t * V;

% Choose acceleration channel
switch lower(accel_channel)
    case 'meas'
        if isfield(R.outputs,'accel_meas') && ~isempty(R.outputs.accel_meas)
            az_s = R.outputs.accel_meas(:);
            accel_used = 'meas';
        else
            warning('No "accel_meas" found. Falling back to "accel_clean".');
            az_s = R.outputs.accel_clean(:);
            accel_used = 'clean';
        end
    case 'clean'
        az_s = R.outputs.accel_clean(:);
        accel_used = 'clean';
    otherwise
        error('Unknown accel_channel="%s"', accel_channel);
end

% Ground truth (for QA later)
x_true  = R.spatial_road(:);
zr_true = R.zr.x(:);

% Spatial sampling
if isfield(R.meta,'fs_spatial')
    fs_spatial = R.meta.fs_spatial;
    dx = R.zr.dx;
else
    dx = mean(diff(x_true),'omitnan');
    fs_spatial = 1/max(dx,eps);
end
nyq_spatial = fs_spatial/2;

%% ======================== INVERSION PARAMETERS ============================
switch lower(param_source)
    case 'from_results'
        P.ms = R.params.ms_corner;
        P.mu = R.params.mu;
        P.ks = R.params.ksc;
        P.cs = R.params.csc;
        P.kt = R.params.ktc;
        P.ct = R.params.ctc;
        P.label = 'from_results (simulation params)';
    
    case 'golden_car'
        P = struct('ms' , 250, ...
                   'mu' , 35, ...
                   'ks' , 63.3e3, ...
                   'cs' , 6.53e3, ...
                   'kt' , 653e3, ...
                   'ct' , 0, ...
                   'label','golden_car (IRI defaults)');

    case 'paper_experimental'
        P = struct('ms' , 400, ...
                   'mu' , 30, ...
                   'ks' , 33700, ...
                   'cs' , 2643.5, ...
                   'kt' , 345400, ...
                   'ct' , 100, ...
                   'label','paper_experimental (Fortuner per-side)');

    case 'paper_numerical'
        P = struct('ms' , 4450, ...
                   'mu' , 420, ...
                   'ks' , 1.0e6, ...
                   'cs' , 1.5e4, ...
                   'kt' , 1.95e6, ...
                   'ct' , 1.0e3, ...
                   'label','paper_numerical (Table 1)');
    otherwise
        error('Unknown param_source="%s"', param_source);
end

ct_nominal = P.ct;
if ct_nominal == 0
    ct_nominal = 100;
end

%% ====== RESONANCES & RECOMMENDED RECONSTRUCTION CUTOFF ====================
f_wh  = (1/(2*pi))*sqrt(P.kt/P.mu);
f_sp  = (1/(2*pi))*sqrt(P.ks/P.ms);

switch lower(lp_rec_rule)
    case 'one_third_wheelhop'
        lp_rec_fc_hz = f_wh / 3;
    otherwise
        error('Unknown lp_rec_rule="%s"', lp_rec_rule);
end
nc_spatial = lp_rec_fc_hz / max(V, eps);

%% ------------------------ QA PLOTS ----------------------------------------
% Figure 1 — Signals
f1 = mkfig('Reconstruction Init — Signals (time domain)', [100 70 1100 700]);
tlo = tiledlayout(f1,3,1,'TileSpacing','compact','Padding','compact');

nexttile(tlo,1);
plot(t, az_s, 'LineWidth', 1.1); grid on; box on; formatAxes();
title(sprintf('Sprung Acceleration a_z^s(t)  (%s)', accel_used));
xlabel('Time (s)'); ylabel('m/s^2'); yline(0,'k:');

% Integrate preview
vzs_raw = cumtrapz(t, az_s);
zss_raw = cumtrapz(t, vzs_raw);
[zss_prev, vzs_prev] = preview_drift_control(zss_raw, vzs_raw, fs, hp_disp_fc_hz, hp_disp_order);

nexttile(tlo,2);
plot(t, vzs_raw, 'Color',[0.25 0.45 0.85]); hold on;
plot(t, vzs_prev, 'r-'); grid on; box on; formatAxes();
title('Integrated Body Velocity — drift preview');
xlabel('Time (s)'); ylabel('m/s'); legend('raw','preview HP');

nexttile(tlo,3);
plot(t, zss_raw, 'Color',[0.25 0.7 0.3]); hold on;
plot(t, zss_prev, 'r-'); grid on; box on; formatAxes();
title('Integrated Body Displacement — drift preview');
xlabel('Time (s)'); ylabel('m'); legend('raw','preview HP');

if fig_export_png, safe_exportgraphics(f1, fullfile(fig_dir,'01_signals_time.png')); end

% Figure 2 — Spectra
f2 = mkfig('Reconstruction Init — Spectra & Cutoffs', [130 90 1100 720]);
tlo2 = tiledlayout(f2,2,1,'TileSpacing','compact','Padding','compact');

nexttile(tlo2,1);
[pp,ff] = welch_psd(az_s, fs);
if ~isempty(pp)
    semilogy(ff, pp, 'LineWidth', 1.2); grid on; box on; formatAxes();
    title('PSD of Body Acceleration  a_z^s(t)'); xlabel('Hz'); ylabel('PSD (m^2/s^4/Hz)');
    add_vline(f_sp,  '--', 'Sprung f_s_p'); add_vline(f_wh,  '--', 'Wheel-hop f_w_h');
    add_vline(lp_rec_fc_hz, '-', sprintf('Rec LP ~ %.2f Hz', lp_rec_fc_hz));
else
    axis off; text(0.1,0.5,'Signal too short for PSD','Units','normalized');
end

nexttile(tlo2,2);
[Zraw,Fraw] = one_sided_fft(zss_raw, fs);
[Zhp, Fhp]  = one_sided_fft(zss_prev, fs);
if ~isempty(Zraw)
    loglog(Fraw, abs(Zraw), 'Color',[0.2 0.6 0.2]); hold on;
    loglog(Fhp,  abs(Zhp) , 'r-'); grid on; box on; formatAxes();
    title('FFT of z_s(t) — raw vs preview HP'); xlabel('Hz'); ylabel('|Z_s(f)| (m)');
else
    axis off; text(0.1,0.5,'Signal too short for FFT','Units','normalized');
end

if fig_export_png
    safe_exportgraphics(f2, fullfile(fig_dir,'02_spectra_and_cutoffs.png'));
end

% Figure 3 — Spatial sampling
f3 = mkfig('Reconstruction Init — Spatial sampling & cutoffs', [160 120 950 360]);
txt = {
    sprintf('V = %.2f m/s (%.1f km/h)', V, V_kmh)
    sprintf('fs = %.1f Hz   |   fs_spatial = %.3f samples/m (Nyquist = %.3f)', fs, fs_spatial, nyq_spatial)
    sprintf('Resonances: f_sp=%.2f Hz, f_wh=%.2f Hz', f_sp, f_wh)
    sprintf('Suggested LP (temporal) = %.2f Hz', lp_rec_fc_hz)
    sprintf('Suggested LP (spatial)  = %.3f cycles/m', nc_spatial)
    };
annotation('textbox',[.06 .22 .88 .7],'String',txt,'FitBoxToText','on',...
    'BackgroundColor',[.97 .97 .98],'FontName','Consolas','FontSize',11);
axis off;
if fig_export_png
    safe_exportgraphics(f3, fullfile(fig_dir,'03_spatial_sampling.png'));
end

%% ===================== SAVE recon configuration ===========================
cfg = struct();
cfg.schema_version = 'recon_cfg_v1';
cfg.fs = fs;
cfg.dt = dt;
cfg.V = V;
cfg.V_kmh = V_kmh;
cfg.dx = dx;
cfg.fs_spatial = fs_spatial;
cfg.nyq_spatial = nyq_spatial;

cfg.hp_disp_fc_hz = hp_disp_fc_hz;
cfg.hp_disp_order = hp_disp_order;

cfg.lp_rec_rule = lp_rec_rule;
cfg.lp_rec_fc_hz = lp_rec_fc_hz;
cfg.lp_rec_order = lp_rec_order;
cfg.nc_spatial = nc_spatial;

cfg.ms = P.ms;
cfg.mu = P.mu;
cfg.ks = P.ks;
cfg.cs = P.cs;
cfg.kt = P.kt;
cfg.ct = ct_nominal;

cfg.f_sp = f_sp;
cfg.f_wh = f_wh;
cfg.param_source = P.label;
cfg.scenario = R.meta.scenario;
cfg.t = t;
cfg.x_t = x_t;
cfg.az_s = az_s;
cfg.accel_used = accel_used;
cfg.accel_unit = 'm/s^2';
cfg.az_clean = R.outputs.accel_clean(:);
cfg.ground.x = x_true;
cfg.ground.zr = zr_true;
cfg.paths.sim_mat = sim_mat_path;
cfg.paths.root = save_folder;
cfg.paths.fig_dir = fig_dir;
cfg.paths.cfg_path = cfg_path;

save(cfg_path,'cfg');
fprintf('Saved: %s\n', cfg_path);

fid = fopen(diag_txt,'w');
if fid>0
    fprintf(fid, 'Reconstruction Init Diagnostics\n');
    fprintf(fid, '--------------------------------\n');
    fprintf(fid, 'Param source : %s\n', cfg.param_source);
    fprintf(fid, 'Accel used   : %s (%s)\n', cfg.accel_used, cfg.accel_unit);
    fprintf(fid, 'fs (Hz)      : %.1f\n', cfg.fs);
    fprintf(fid, 'V (m/s)      : %.3f   (%.1f km/h)\n', cfg.V, cfg.V_kmh);
    fprintf(fid, 'Resonances   : f_sp=%.2f Hz, f_wh=%.2f Hz\n', cfg.f_sp, cfg.f_wh);
    fprintf(fid, 'LP cutoff    : temporal=%.2f Hz, spatial=%.3f cyc/m\n', cfg.lp_rec_fc_hz, cfg.nc_spatial);
    fprintf(fid, 'HP disp prev : fc=%.2f Hz, order=%d\n', cfg.hp_disp_fc_hz, cfg.hp_disp_order);
    fprintf(fid, 'Spatial fs   : %.3f samples/m (Nyquist=%.3f)\n', cfg.fs_spatial, cfg.nyq_spatial);
    fclose(fid);
end
fprintf('Saved: %s\n', diag_txt);

%% =============================== FUNCTIONS ===============================
function f = mkfig(name, pos)
    f = figure('Name', name, 'Color', 'w', 'Position', pos);
end

function formatAxes(ax)
    if nargin < 1
        ax = gca;
    end
    ax.LineWidth = 1.0;
    ax.FontName = 'Calibri';
    ax.FontSize = 11;
    grid(ax,'on'); box(ax,'on');
end

function add_vline(x, style, labelStr)
    if isempty(x) || ~isfinite(x)
        return;
    end
    yl = ylim;
    hold on;
    plot([x x], yl, style, 'Color',[.15 .15 .15], 'LineWidth', 1.0);
    text(x, yl(2), sprintf('  %s (%.2f)', labelStr, x), ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'Rotation',90,'FontSize',9,'Color',[.15 .15 .15]);
end

function [Pxx,F] = welch_psd(x, fs)
    % Robust Welch PSD: nperseg ≤ N, avoids short-signal errors.
    x = x(:) - mean(x(:),'omitnan');
    N = numel(x);
    if N < 8
        Pxx = []; F = [];
        return;
    end
    nper = min(N, max(128, 2^floor(log2(max(32, floor(N/4))))));
    win  = hamming(nper);
    nover= floor(nper/2);
    nfft = 2^nextpow2(max(nper, min(8192, N)));
    [Pxx,F] = pwelch(x, win, nover, nfft, fs, 'onesided');
end

function [Xmag, F] = one_sided_fft(x, fs)
    x = x(:) - mean(x(:),'omitnan');
    N = numel(x);
    if N < 2
        Xmag = []; F = [];
        return;
    end
    Nfft = max(2, 2^nextpow2(N));
    X = fft(x, Nfft);
    X = X(1:floor(Nfft/2)+1);
    Xmag = X / N;
    F = fs*(0:floor(Nfft/2))/Nfft;
end

function [z_hp, vz_hp] = preview_drift_control(z_raw, vz_raw, fs, fc, order)
    if fc<=0
        z_hp  = z_raw;
        vz_hp = vz_raw;
        return;
    end
    Wn = max(min(fc/(fs/2), 0.999), 1e-6);
    [b,a] = butter(order, Wn, 'high');
    z_hp  = filtfilt(b,a,z_raw);
    vz_hp = filtfilt(b,a,vz_raw);
end

function safe_exportgraphics(figHandle, outPath)
    try
        exportgraphics(figHandle, outPath, 'Resolution', 180);
    catch
        % Fallback for older MATLAB
        [folder, base, ~] = fileparts(outPath);
        if ~exist(folder,'dir')
            mkdir(folder);
        end
        saveas(figHandle, fullfile(folder, [base '.png']));
    end
end