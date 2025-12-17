%% Load Real Data and Preprocessing
% -------------------------------------------------------------------------
% PURPOSE:
%   Imports pre-converted MAT data, selects a specific signal, and applies
%   preprocessing STRICTLY according to the reference thesis:
%     1. Zero-Mean (DC Removal)
%     2. Wavelet Denoising (Coif5, Level 5, Soft Threshold)
%
%   UPDATE: Added toggle to enable/disable cleaning steps.
%   UPDATE: Fixed unit consistency in plots, improved error reporting,
%           separated accel_meas vs accel_clean, and made PSD more robust.
%
% OUTPUT:
%   ./00_Outputs/01_Simulation/SimulationData/simulation_results.mat
% -------------------------------------------------------------------------

close all; clear; clc;

%% ======================= 0) USER HYPERPARAMETERS ========================
% --- Data Source (MAT File) ---
P.input_mat_file  = fullfile('Real_Data', 'Acceleration_Data.mat');

% --- Selection ---
P.target_sheet    = 'Sheet1';
P.target_column   = 'sig1_new';
% P.target_column   = 'sig2_shifted';

% --- Experiment Parameters ---
P.fs              = 200;                 % Sampling frequency (Hz)
P.speed_kmh       = 30;                  % Average vehicle speed (km/h)
P.units_are_g     = true;                % Set TRUE if input is in 'g'

% --- Processing Control ---
% Set TRUE to apply Zero-Mean and Wavelet Denoising (Thesis Method).
% Set FALSE to use raw data (only unit conversion will be applied).
P.enable_cleaning = false;                

% --- Vehicle Assumption ---
P.veh_preset      = 'Golden_Car';

% --- Preprocessing Settings (Thesis Section 4.2) ---
% Thesis Table 4-3: Smartphone data uses 'coif5', level 5, Soft threshold.
P.wavelet_name    = 'coif5';
P.wavelet_level   = 5;
P.threshold_mode  = 's';                 % 's' for soft, 'h' for hard

% --- Plotting ---
P.fig_size        = [100, 100, 1000, 800];
P.save_plots      = true;

%% 1) Load Data
fprintf('Loading data from: %s ...\n', P.input_mat_file);
if ~isfile(P.input_mat_file)
    error('File not found. Run "ConvertExcelToMat.m" first.');
end

loadedData = load(P.input_mat_file);
fn = fieldnames(loadedData);
AllSignals = loadedData.(fn{1});

if ~isfield(AllSignals, P.target_sheet)
    error('Sheet "%s" not found.', P.target_sheet);
end
targetTable = AllSignals.(P.target_sheet);

vars = targetTable.Properties.VariableNames;
if ~ismember(P.target_column, vars)
    error('Column "%s" not found in %s.', P.target_column, P.target_sheet);
end

raw_accel = targetTable.(P.target_column);
fprintf('Loaded: %s -> %s (%d samples)\n', P.target_sheet, P.target_column, length(raw_accel));

N = length(raw_accel);
dt = 1/P.fs;
raw_time = (0:N-1)' * dt;

%% 2) Pre-Processing
accel_proc = raw_accel(:);

% --- A. Unit Conversion (Always run to ensure m/s^2) ---
if P.units_are_g
    accel_proc    = accel_proc * 9.80665;
    raw_accel_ms2 = raw_accel * 9.80665;   % for plotting/PSD
    fprintf('Converted units to m/s^2.\n');
else
    raw_accel_ms2 = raw_accel;            % already in m/s^2
end

% --- B. Cleaning Steps (Controlled by P.enable_cleaning) ---
if P.enable_cleaning
    fprintf('Preprocessing ENABLED (Thesis Method)...\n');
    
    % 1. Zero-Mean (Thesis Section 3.5.3, Eq 3-6)
    signal_mean   = mean(accel_proc);
    accel_zeromean = accel_proc - signal_mean;
    fprintf('  - Applied Zero-Mean (Bias: %.4f)\n', signal_mean);

    % 2. Wavelet Denoising (Thesis Section 3.5.4 & 4.2.3)
    fprintf('  - Applying Denoising: %s, Level %d, %s-Threshold...\n', ...
        P.wavelet_name, P.wavelet_level, P.threshold_mode);

    try
        % Decompose signal using DWT
        [C, L] = wavedec(accel_zeromean, P.wavelet_level, P.wavelet_name);

        % Calculate adaptive threshold (Universal Threshold / VisuShrink)
        cD1   = detcoef(C, L, 1);
        sigma = median(abs(cD1)) / 0.6745;
        thr   = sigma * sqrt(2 * log(N));

        % Perform thresholding on detail coefficients
        % 'keepapp' = 1 ensures we do NOT threshold the approximation (Low Freq)
        [accel_denoised, ~, ~, ~] = wdencmp('gbl', C, L, P.wavelet_name, ...
            P.wavelet_level, thr, P.threshold_mode, 1);

        accel_final = accel_denoised;
        fprintf('  - Denoising complete.\n');

    catch
        warning('Wavelet Toolbox error: %s');
        fprintf('  - Falling back to Zero-Mean signal only.\n');
        accel_final = accel_zeromean;
    end
else
    fprintf('Preprocessing DISABLED. Using Raw Data (Unit Converted).\n');
    % If cleaning is disabled, pass the raw (but unit-converted) signal
    accel_final    = accel_proc;
    accel_zeromean = accel_proc; 
end

% --- B. Cleaning Steps (Controlled by P.enable_cleaning) ---
if P.enable_cleaning
    fprintf('Preprocessing ENABLED (Thesis Method)...\n');
    
    % 1. Zero-Mean (Thesis Section 3.5.3, Eq 3-6)
    signal_mean    = mean(accel_proc);
    accel_zeromean = accel_proc - signal_mean;
    fprintf('  - Applied Zero-Mean (Bias: %.4f)\n', signal_mean);

    % 2. Wavelet Denoising (Thesis Section 3.5.4 & 4.2.3)
    fprintf('  - Applying Denoising: %s, Level %d, %s-Threshold...\n', ...
        P.wavelet_name, P.wavelet_level, P.threshold_mode);

    try
        % Decompose signal using DWT
        [C, L] = wavedec(accel_zeromean, P.wavelet_level, P.wavelet_name);

        % Calculate adaptive threshold (Universal Threshold / VisuShrink)
        cD1   = detcoef(C, L, 1);
        sigma = median(abs(cD1)) / 0.6745;
        thr   = sigma * sqrt(2 * log(N));

        % Perform thresholding on detail coefficients
        % 'keepapp' = 1 ensures we do NOT threshold the approximation (Low Freq)
        [accel_denoised, ~, ~, ~] = wdencmp('gbl', C, L, P.wavelet_name, ...
            P.wavelet_level, thr, P.threshold_mode, 1);

        accel_final = accel_denoised;
        fprintf('  - Denoising complete.\n');

    catch
        warning('Wavelet Toolbox error: %s');
        fprintf('  - Falling back to Zero-Mean signal only.\n');
        accel_final = accel_zeromean;
    end
else
    fprintf('Preprocessing DISABLED. Using Raw Data (Unit Converted).\n');
    % If cleaning is disabled, pass the raw (but unit-converted) signal
    accel_final    = accel_proc;
    accel_zeromean = accel_proc; 
end

% Residual (noise removed by denoising)
noise_residual = accel_zeromean - accel_final;

% SNR estimate: treat accel_final as "signal" and noise_residual as "noise"
signal_power = mean(accel_final.^2);
noise_power  = mean(noise_residual.^2);
SNR_dB       = 10*log10(signal_power / noise_power);
fprintf('Estimated SNR (denoised vs residual): %.2f dB\n', SNR_dB);

%% 3) Package Results
results = struct();
results.meta.vehicle_model      = 'Real_World_Data';
results.meta.scenario           = ['Recorded_' P.target_column];
results.meta.road_class         = 'Unknown (Blind)';
results.meta.V_kmh              = P.speed_kmh;
results.meta.fs                 = P.fs;
results.meta.dt                 = dt;
results.meta.timestamp          = char(datetime('now'));
results.meta.source_file        = P.input_mat_file;
results.meta.cleaning_applied   = P.enable_cleaning;
results.meta.wavelet_name       = P.wavelet_name;
results.meta.wavelet_level      = P.wavelet_level;
results.meta.threshold_mode     = P.threshold_mode;

% Vehicle Params
results.params = get_vehicle_dummy_params(P.veh_preset);

% Vectors
V_ms          = P.speed_kmh / 3.6;
results.time  = raw_time(:);
results.x_m   = raw_time(:) * V_ms;

% Outputs
results.outputs.accel_meas  = accel_proc;   % unit-converted measurement
results.outputs.accel_clean = accel_final;  % final processed signal
results.outputs.heave_m     = zeros(N, 1);

% Dummy Ground Truth
results.spatial_road    = results.x_m;
results.zr.x            = zeros(N, 1);
results.zr.base_x       = zeros(N, 1);
results.zr.t            = zeros(N, 1);
results.zr.dz_dt        = zeros(N, 1);
results.zr.dt           = dt;
results.zr.dx           = V_ms * dt;
results.meta.fs_spatial = 1 / results.zr.dx;

% Fillers
results.inputs.U      = zeros(N, 2);
results.inputs.labels = {'zr','dzr'};
results.states.X      = zeros(N, 4);
results.iso8608       = [];
results.smartphone    = [];
% Metrics
results.metrics.SNR_dB = SNR_dB;

%% 4) Save
out_dir = fullfile('00_Outputs', '01_Simulation', 'SimulationData');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save_path = fullfile(out_dir, 'simulation_results.mat');
save(save_path, 'results', '-v7.3');
fprintf('\nSUCCESS! Data saved to:\n  %s\n', save_path);

%% 5) Diagnostics (Thesis-style Validation)
nrows = 4;  % time-domain (2), residual (1), PSD (1)
fig = figure('Name', ['Preprocessing: ' P.target_column], ...
             'Color', 'w', 'Position', P.fig_size);
tiledlayout(nrows, 1, 'Padding', 'compact');

% --- Plot 1: Raw vs. Zero-Mean (or Final if cleaning off) ---
nexttile;
plot(raw_time, raw_accel_ms2, 'Color', [0.6 0.6 0.6], 'DisplayName', 'Raw Input');
hold on;
if P.enable_cleaning
    plot(raw_time, accel_zeromean, 'b', 'LineWidth', 0.8, 'DisplayName', 'Zero-Mean');
    title('Step 1: Zero-Mean (DC Removal)');
else
    plot(raw_time, accel_zeromean, 'b', 'LineWidth', 0.8, 'DisplayName', 'Unit-Converted');
    title('Step 1: Raw Data (Cleaning Disabled)');
end
legend; grid on; xlim([0, max(raw_time)]); ylabel('Accel (m/s^2)');

% --- Plot 2: Zero-Mean vs. Denoised (or Final) ---
nexttile;
if P.enable_cleaning
    plot(raw_time, accel_zeromean, 'Color', [0.6 0.6 1.0], ...
        'DisplayName', 'Noisy (Zero-Mean)');
    hold on;
    plot(raw_time, accel_final, 'r', 'LineWidth', 1.2, ...
        'DisplayName', 'Denoised (Coif5)');
    title(['Step 2: Wavelet Denoising (' P.wavelet_name ...
           ' L' num2str(P.wavelet_level) ')']);
else
    plot(raw_time, accel_final, 'k', 'LineWidth', 1.0, ...
        'DisplayName', 'Unit-Converted');
    title('Step 2: Skipped (Cleaning Disabled)');
end
legend; grid on; xlim([0, max(raw_time)]); ylabel('Accel (m/s^2)');

% --- Plot 3: Residual (Zero-Mean - Denoised) + SNR annotation ---
nexttile;
if P.enable_cleaning
    plot(raw_time, noise_residual, 'm', 'LineWidth', 1.0, ...
        'DisplayName', 'Residual (Zero-Mean - Denoised)');
    title(sprintf('Step 3: Removed Noise (Residual)  |  SNR: approx %.2f dB', SNR_dB));
    ylabel('Accel (m/s^2)');
    grid on; xlim([0, max(raw_time)]);
    legend;
else
    axis off;
    text(0.5, 0.5, 'Residual plot skipped (cleaning disabled)', ...
        'HorizontalAlignment', 'center', 'FontSize', 11);
end

% --- Plot 4: PSD Comparison (Raw, Zero-Mean, Denoised) ---
nexttile;
% Make PSD robust for short signals
nfft     = min(1024, N);
win      = hamming(min(1024, N));
noverlap = floor(length(win)/2);

[pxx_raw,   f] = pwelch(raw_accel_ms2, win, noverlap, nfft, P.fs);
[pxx_zero, ~]  = pwelch(accel_zeromean, win, noverlap, nfft, P.fs);
[pxx_final, ~] = pwelch(accel_final,   win, noverlap, nfft, P.fs);

loglog(f, pxx_raw,  'Color', [0.3 0.3 1.0], 'DisplayName', 'Raw Input');
hold on;
loglog(f, pxx_zero, 'Color', [0.6 0.6 1.0], 'DisplayName', 'Zero-Mean');
loglog(f, pxx_final,'r', 'LineWidth', 1.5, 'DisplayName', 'Output Signal');
title('Frequency Domain Check');
xlabel('Hz'); ylabel('PSD');
legend; grid on; xlim([0.5, 100]);

if P.save_plots
    exportgraphics(fig, fullfile(out_dir, 'preprocessing_check.png'), 'Resolution', 300);
end

%% HELPER
function params = get_vehicle_dummy_params(preset_name)
    switch lower(preset_name)
        case 'sport_sedan'
            ms = 240; mu = 40; ks = 85e3; cs = 9.0e3; kt = 800e3; ct = 200;
        case 'suv'
            ms = 320; mu = 45; ks = 65e3; cs = 7.5e3; kt = 600e3; ct = 150;
        case 'compact'
            ms = 300; mu = 40; ks = 75e3; cs = 8.0e3; kt = 720e3; ct = 220;
        otherwise  % 'golden_car' or anything else
            ms = 250; mu = 35; ks = 63.3e3; cs = 6.53e3; kt = 653e3; ct = 0;
    end
    params.ms_corner = ms; 
    params.mu        = mu;
    params.ksc       = ks; 
    params.csc       = cs;
    params.ktc       = kt; 
    params.ctc       = ct;
end
