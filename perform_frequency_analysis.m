
%% 1-b) Perform Comprehensive Frequency Analysis of Vehicle Acceleration
% -------------------------------------------------------------------------
% PURPOSE:
% This script analyzes the frequency content of the simulated vehicle's
% vertical acceleration. It uses techniques like the Fast Fourier Transform
% (FFT), Power Spectral Density (PSD), and spectrograms to identify
% dominant frequencies in the vehicle's response and compares them to the
% theoretical natural frequencies of the vehicle model.
%
% CORE LOGIC (How it works):
% 1.  Loads the 'simulation_results.mat' file to get the time-history
%     data, vehicle parameters, and simulation settings.
%
% 2.  Calculates the theoretical natural frequencies (body and wheel-hop
%     modes) using the exact vehicle parameters loaded from the simulation
%     file. This ensures the theoretical analysis matches the simulation.
%
% 3.  Performs an FFT on the acceleration signal to find its raw magnitude
%     and phase spectrum, identifying the most dominant frequency components.
%
% 4.  Calculates the PSD using Welch's method, which provides a smoothed,
%     more robust estimate of how the signal's power is distributed
%     across different frequencies.
%
% 5.  Generates a spectrogram, which visualizes how the frequency content
%     of the acceleration signal changes over the duration of the simulation,
%     making it possible to link frequency events to specific moments in time.
%
% 6.  Creates a series of detailed plots summarizing the time-domain stats,
%     FFT, PSD, and spectrogram results, saving them as images. It also
%     saves all calculated analysis data into a new .mat file and a summary
%     .txt report.
%
% INPUTS:
% - '01_Quarter Car Modeling/simulation_results.mat': The output file from
%   the main vehicle simulation.
%
% OUTPUTS:
% - '01_2_Frequency_Analysis/frequency_analysis_results.mat': A file
%   containing all calculated spectral data and metrics.
% - A .txt summary report of the analysis.
% - Multiple .png figures visualizing the results.
% -------------------------------------------------------------------------

close all; clc;

%% --------------- Load simulation results (consistent with 01_Simulation) ---------------
sim_file = fullfile('00_Outputs','01_Simulation','SimulationData','simulation_results.mat');
assert(exist(sim_file,'file')==2, 'Missing results MAT: %s. Run the simulation script first.', sim_file);

S = load(sim_file); results = S.results;

% Canonical variables
t_sol   = results.time(:);
V       = results.meta.V_kmh / 3.6;          % m/s
Vehicle_Speed = results.meta.V_kmh;          % km/h (display)
Sampling_Freq = results.meta.fs;

spatial_road  = results.spatial_road(:);
z_r_x         = results.zr.x(:);             % spatial road
z_r_t         = results.zr.t(:);             % time-domain road (not strictly needed here)

Y_sol               = results.states.X;      % [ys dys yus dyus]
Sprung_Mass_Disp    = results.outputs.heave_m(:);
Unsprung_Mass_Disp  = Y_sol(:,3);

% Choose noisy (measured) accel for analysis; clean is available too
sprung_mass_accel      = results.outputs.accel_meas(:);
sprung_mass_accel_raw  = results.outputs.accel_clean(:);

road_class                = results.meta.road_class;
damage_scenario_selection = results.meta.scenario;
if isfield(results,'meta') && isfield(results.meta,'solver_method')
    solver_method = results.meta.solver_method;
else
    solver_method = 'unknown';
end

%% --------------- Analysis parameters ---------------
% Frequency analysis parameters
max_freq_display = 100;                  % Hz - maximum frequency to display
min_freq_display = 0.1;                 % Hz - minimum frequency to display
window_length = 1.0;                    % seconds - window length for PSD
overlap_ratio = 0.75;                   % 75% overlap between windows

% Vehicle parameters for theoretical analysis (from results)
ms = results.params.ms_corner;
mus = results.params.mu;
ks = results.params.ksc;
cs = results.params.csc;
kt = results.params.ktc;
ct = results.params.ctc;

% Natural frequencies & damping from quarter-car state matrix
A_chk = [0 1 0 0;
         -ks/ms -cs/ms  ks/ms   cs/ms;
          0 0 0 1;
          ks/mus cs/mus -(ks+kt)/mus -(cs+ct)/mus];
ev = eig(A_chk);
wn_all   = abs(imag(ev));
zeta_all = -real(ev)./sqrt(real(ev).^2 + imag(ev).^2 + eps);
f_all    = wn_all/(2*pi);

[fsorted, idxm] = sort(f_all);
[~, ia] = unique(round(fsorted,2),'stable');
fu = fsorted(ia); zu = zeta_all(idxm(ia));
assert(numel(fu)>=2,'Could not identify two modes from parameters.');

fn1  = fu(1);   % lower = body mode
fn2  = fu(2);   % higher = wheel-hop
zeta1 = zu(1);
zeta2 = zu(2);

%% --------------- Create output folder ---------------
output_folder = fullfile('00_Outputs', '01_Simulation', 'FrequencyAnalysis');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% --------------- Data preparation ---------------
% Ensure all data is column vectors
t_sol = t_sol(:);
sprung_mass_accel = sprung_mass_accel(:);
Sprung_Mass_Disp = Sprung_Mass_Disp(:);
Unsprung_Mass_Disp = Unsprung_Mass_Disp(:);

% Calculate derived signals
dt = mean(diff(t_sol));
sprung_mass_vel = gradient(Sprung_Mass_Disp, dt);
unsprung_mass_vel = gradient(Unsprung_Mass_Disp, dt);

% Calculate road input velocity
road_vel = gradient(z_r_x, spatial_road) * V;

% Signal statistics
accel_mean = mean(sprung_mass_accel);
accel_std = std(sprung_mass_accel);
accel_rms = rms(sprung_mass_accel);
accel_max = max(abs(sprung_mass_accel));

fprintf('Acceleration Signal Statistics:\n');
fprintf('  Mean: %.4f m/s²\n', accel_mean);
fprintf('  Std: %.4f m/s²\n', accel_std);
fprintf('  RMS: %.4f m/s²\n', accel_rms);
fprintf('  Max: %.4f m/s²\n', accel_max);
fprintf('  Duration: %.2f s\n', t_sol(end) - t_sol(1));
fprintf('  Samples: %d\n', length(sprung_mass_accel));
fprintf('\n');

%% --------------- Frequency Analysis ---------------
fprintf('Performing frequency analysis...\n');

% FFT Analysis
N = length(sprung_mass_accel);
N_fft = 2^nextpow2(N);
freq_vector = (0:N_fft/2-1) * Sampling_Freq / N_fft;
fft_result = fft(sprung_mass_accel, N_fft);
fft_magnitude = abs(fft_result(1:N_fft/2));
fft_phase = angle(fft_result(1:N_fft/2));

% Power Spectral Density (PSD)
window_samples = round(window_length * Sampling_Freq);
overlap_samples = round(overlap_ratio * window_samples);
window_samples = min(length(sprung_mass_accel), window_samples);
overlap_samples = max(0, min(overlap_samples, window_samples-1));
[psd, f_psd] = pwelch(sprung_mass_accel, hann(window_samples,'periodic'), overlap_samples, N_fft, Sampling_Freq);

% Spectrogram
[~, F, T, P] = spectrogram(sprung_mass_accel, hann(window_samples,'periodic'), overlap_samples, N_fft, Sampling_Freq, 'yaxis');

% Filter frequency range for display
freq_idx = (freq_vector >= min_freq_display) & (freq_vector <= max_freq_display);
freq_display = freq_vector(freq_idx);
fft_magnitude_display = fft_magnitude(freq_idx);
fft_phase_display = fft_phase(freq_idx);

f_psd_idx = (f_psd >= min_freq_display) & (f_psd <= max_freq_display);
f_psd_display = f_psd(f_psd_idx);
psd_display = psd(f_psd_idx);

F_idx = (F >= min_freq_display) & (F <= max_freq_display);
F_display = F(F_idx);
P_display = P(F_idx, :);

% Find peak frequencies
[peaks, peak_locs] = findpeaks(fft_magnitude_display, 'MinPeakHeight', max(fft_magnitude_display)*0.1, 'MinPeakDistance', 5);
peak_frequencies = freq_display(peak_locs);
peak_magnitudes = peaks;

% Find dominant frequencies
[~, dominant_idx] = max(fft_magnitude_display);
dominant_frequency = freq_display(dominant_idx);

fprintf('Frequency Analysis Results:\n');
fprintf('  Dominant frequency: %.2f Hz\n', dominant_frequency);
fprintf('  Theoretical body mode: %.2f Hz\n', fn1);
fprintf('  Theoretical wheel hop: %.2f Hz\n', fn2);
fprintf('  Number of peaks detected: %d\n', length(peak_frequencies));
fprintf('\n');

%% --------------- Figure 1: Time Domain Analysis ---------------
fig1 = figure('Name', 'Time Domain Analysis', 'Color', 'w', 'Position', [50, 50, 1200, 800]);

% Create tiled layout
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Acceleration signal
ax1 = nexttile(1, [1, 2]);
hold(ax1, 'on');
grid(ax1, 'on'); grid(ax1, 'minor');

plot(ax1, t_sol, sprung_mass_accel, 'b-', 'LineWidth', 1.5);
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Acceleration (m/s²)');
title(ax1, sprintf('Acceleration Signal - %s (Class %s)', strrep(damage_scenario_selection,'_',' '), road_class));
xlim(ax1, [t_sol(1), t_sol(end)]);

% Add statistics text
stats_text = sprintf('Mean: %.3f m/s²\nStd: %.3f m/s²\nRMS: %.3f m/s²\nMax: %.3f m/s²', ...
                     accel_mean, accel_std, accel_rms, accel_max);
text(ax1, 0.02, 0.98, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'k');

% Subplot 2: Acceleration histogram
ax2 = nexttile(3);
hold(ax2, 'on');
grid(ax2, 'on'); grid(ax2, 'minor');

histogram(ax2, sprung_mass_accel, 50, 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha', 0.7);
xlabel(ax2, 'Acceleration (m/s²)');
ylabel(ax2, 'Frequency');
title(ax2, 'Acceleration Distribution');

% Add normal distribution overlay
x_norm = linspace(min(sprung_mass_accel), max(sprung_mass_accel), 100);
y_norm = normpdf(x_norm, mean(sprung_mass_accel), std(sprung_mass_accel));
y_norm = y_norm * max(histcounts(sprung_mass_accel, 50)) / max(y_norm);
plot(ax2, x_norm, y_norm, 'r-', 'LineWidth', 2);

% Subplot 3: Cumulative RMS
ax3 = nexttile(4);
hold(ax3, 'on');
grid(ax3, 'on'); grid(ax3, 'minor');

cumulative_rms = sqrt( cumtrapz(t_sol, sprung_mass_accel.^2) ./ max(t_sol, eps) );
plot(ax3, t_sol, cumulative_rms, 'r-', 'LineWidth', 1.5);
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'RMS Acceleration (m/s²)');
title(ax3, 'Cumulative RMS Acceleration');

% Save Figure 1
saveas(fig1, fullfile(output_folder, '01_time_domain_analysis.png'));

%% --------------- Figure 2: Frequency Domain Analysis ---------------
fig2 = figure('Name', 'Frequency Domain Analysis', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

% Create tiled layout
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: FFT Magnitude
ax1 = nexttile(1);
hold(ax1, 'on');
grid(ax1, 'on'); grid(ax1, 'minor');

plot(ax1, freq_display, fft_magnitude_display, 'b-', 'LineWidth', 2);

% Add theoretical natural frequencies
xline(ax1, fn1, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Body Mode: %.1f Hz', fn1));
xline(ax1, fn2, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Wheel Hop: %.1f Hz', fn2));

% Mark peak frequencies
if ~isempty(peak_frequencies)
    plot(ax1, peak_frequencies, peak_magnitudes, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

xlabel(ax1, 'Frequency (Hz)');
ylabel(ax1, 'Magnitude');
title(ax1, 'FFT Magnitude');
legend(ax1, {'FFT', 'Body Mode', 'Wheel Hop', 'Peaks'}, 'Location', 'best');
xlim(ax1, [min_freq_display, max_freq_display]);

% Subplot 2: Power Spectral Density
ax2 = nexttile(2);
hold(ax2, 'on');
grid(ax2, 'on'); grid(ax2, 'minor');

semilogy(ax2, f_psd_display, psd_display, 'b-', 'LineWidth', 2);

% Add theoretical natural frequencies
xline(ax2, fn1, 'r--', 'LineWidth', 2);
xline(ax2, fn2, 'g--', 'LineWidth', 2);

xlabel(ax2, 'Frequency (Hz)');
ylabel(ax2, 'PSD (m²/s⁴/Hz)');
title(ax2, 'Power Spectral Density');
xlim(ax2, [min_freq_display, max_freq_display]);

% Subplot 3: FFT Phase
ax3 = nexttile(3);
hold(ax3, 'on');
grid(ax3, 'on'); grid(ax3, 'minor');

plot(ax3, freq_display, fft_phase_display*180/pi, 'b-', 'LineWidth', 1.5);

xlabel(ax3, 'Frequency (Hz)');
ylabel(ax3, 'Phase (degrees)');
title(ax3, 'FFT Phase');
xlim(ax3, [min_freq_display, max_freq_display]);
ylim(ax3, [-180, 180]);

% Subplot 4: Peak Analysis
ax4 = nexttile(4);
hold(ax4, 'on');
grid(ax4, 'on'); grid(ax4, 'minor');

if ~isempty(peak_frequencies)
    bar(ax4, peak_frequencies, peak_magnitudes, 'FaceColor', 'blue', 'EdgeColor', 'black');
    
    % Add theoretical frequencies
    xline(ax4, fn1, 'r--', 'LineWidth', 2);
    xline(ax4, fn2, 'g--', 'LineWidth', 2);
    
    xlabel(ax4, 'Frequency (Hz)');
    ylabel(ax4, 'Magnitude');
    title(ax4, 'Peak Frequency Analysis');
    legend(ax4, {'Peaks', 'Body Mode', 'Wheel Hop'}, 'Location', 'best');
else
    text(ax4, 0.5, 0.5, 'No significant peaks detected', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
    title(ax4, 'Peak Frequency Analysis');
end

% Save Figure 2
saveas(fig2, fullfile(output_folder, '02_frequency_domain_analysis.png'));

%% --------------- Figure 3: Spectrogram Analysis ---------------
fig3 = figure('Name', 'Spectrogram Analysis', 'Color', 'w', 'Position', [150, 150, 1200, 600]);

% Create tiled layout
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Spectrogram
ax1 = nexttile(1);
hold(ax1, 'on');

imagesc(ax1, T, F_display, 10*log10(P_display + eps));
colorbar(ax1);
colormap(ax1, 'jet');

% Add theoretical natural frequencies
yline(ax1, fn1, 'r--', 'LineWidth', 2);
yline(ax1, fn2, 'g--', 'LineWidth', 2);

xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Frequency (Hz)');
title(ax1, 'Spectrogram (dB)');
ylim(ax1, [min_freq_display, max_freq_display]);

% Subplot 2: Cumulative Power
ax2 = nexttile(2);
hold(ax2, 'on');
grid(ax2, 'on'); grid(ax2, 'minor');

% Calculate cumulative power
cumulative_power = cumsum(psd_display);
cumulative_power_norm = cumulative_power / max(cumulative_power) * 100;

plot(ax2, f_psd_display, cumulative_power_norm, 'b-', 'LineWidth', 2);

xlabel(ax2, 'Frequency (Hz)');
ylabel(ax2, 'Cumulative Power (%)');
title(ax2, 'Cumulative Power Distribution');
xlim(ax2, [min_freq_display, max_freq_display]);
ylim(ax2, [0, 100]);

% Save Figure 3
saveas(fig3, fullfile(output_folder, '03_spectrogram_analysis.png'));

%% --------------- Figure 4: Theoretical Analysis ---------------
fig4 = figure('Name', 'Theoretical Analysis', 'Color', 'w', 'Position', [200, 200, 1200, 800]);

% Create tiled layout
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Frequency Response Comparison
ax1 = nexttile(1);
hold(ax1, 'on');
grid(ax1, 'on'); grid(ax1, 'minor');

% Plot theoretical frequency response
freq_theory = linspace(min_freq_display, max_freq_display, 1000);
s = 1i * 2 * pi * freq_theory;

% Full quarter-car FRF: road elevation zr -> body vertical acceleration \ddot{z}_b
A = [0 1 0 0;
     -ks/ms -cs/ms  ks/ms   cs/ms;
      0 0 0 1;
      ks/mus cs/mus -(ks+kt)/mus -(cs+ct)/mus];
B = [0 0; 0 0; 0 0; kt/mus ct/mus];
Cacc = [-ks/ms, -cs/ms,  ks/ms,  cs/ms];  % \ddot{z}_b = Cacc*x
Dacc = [0 0];

Hmag = zeros(size(freq_theory));
for ii = 1:numel(freq_theory)
    jw = 1i*2*pi*freq_theory(ii);
    Hmag(ii) = abs( Cacc / (jw*eye(4) - A) * B(:,1) ); % zr -> \ddot{z}_b
end

plot(ax1, freq_theory, Hmag, 'm-', 'LineWidth', 2, 'DisplayName', 'FRF | zr \rightarrow \ddot{z}_b');

xlabel(ax1, 'Frequency (Hz)');
ylabel(ax1, 'Magnitude');
title(ax1, 'Theoretical Frequency Response');
legend(ax1, 'Location', 'best');
xlim(ax1, [min_freq_display, max_freq_display]);

% Subplot 2: Modal Analysis
ax2 = nexttile(2);
hold(ax2, 'on');
grid(ax2, 'on'); grid(ax2, 'minor');

% Plot theoretical mode shapes
freq_theory = linspace(0, max_freq_display, 1000);
body_mode = 1./(1 - (freq_theory/fn1).^2 + 1i*2*zeta1*(freq_theory/fn1));
wheel_mode = 1./(1 - (freq_theory/fn2).^2 + 1i*2*zeta2*(freq_theory/fn2));

plot(ax2, freq_theory, abs(body_mode), 'r-', 'LineWidth', 2, 'DisplayName', 'Body Mode');
plot(ax2, freq_theory, abs(wheel_mode), 'g-', 'LineWidth', 2, 'DisplayName', 'Wheel Mode');

xlabel(ax2, 'Frequency (Hz)');
ylabel(ax2, 'Magnitude');
title(ax2, 'Theoretical Mode Shapes');
legend(ax2, 'Location', 'best');
xlim(ax2, [0, max_freq_display]);

% Subplot 3: Comparison with Measured Data
ax3 = nexttile(3);
hold(ax3, 'on');
grid(ax3, 'on'); grid(ax3, 'minor');

% Normalize both for comparison
fft_norm = fft_magnitude_display / max(fft_magnitude_display);
body_mode_norm = abs(body_mode) / max(abs(body_mode));

plot(ax3, freq_display, fft_norm, 'b-', 'LineWidth', 2, 'DisplayName', 'Measured FFT');
plot(ax3, freq_theory, body_mode_norm, 'r--', 'LineWidth', 2, 'DisplayName', 'Theoretical Body Mode');

xlabel(ax3, 'Frequency (Hz)');
ylabel(ax3, 'Normalized Magnitude');
title(ax3, 'Measured vs Theoretical');
legend(ax3, 'Location', 'best');
xlim(ax3, [min_freq_display, max_freq_display]);

% Subplot 4: Analysis Summary
ax4 = nexttile(4);
hold(ax4, 'on');
grid(ax4, 'on'); grid(ax4, 'minor');

% Create analysis summary table
summary_text = {
    'FREQUENCY ANALYSIS SUMMARY';
    '========================';
    '';
    sprintf('Dominant Frequency: %.2f Hz', dominant_frequency);
    sprintf('Theoretical Body Mode: %.2f Hz', fn1);
    sprintf('Theoretical Wheel Hop: %.2f Hz', fn2);
    '';
    sprintf('Body Mode Error: %.2f%%', abs(dominant_frequency - fn1)/fn1*100);
    sprintf('Wheel Hop Error: %.2f%%', abs(dominant_frequency - fn2)/fn2*100);
    '';
    sprintf('Peak Frequencies:');
};

if ~isempty(peak_frequencies)
    for i = 1:min(3, length(peak_frequencies))
        summary_text{end+1} = sprintf('  %.2f Hz', peak_frequencies(i));
    end
else
    summary_text{end+1} = '  None detected';
end

summary_text{end+1} = '';
summary_text{end+1} = 'ANALYSIS PARAMETERS';
summary_text{end+1} = '==================';
summary_text{end+1} = sprintf('Sampling Rate: %.0f Hz', Sampling_Freq);
summary_text{end+1} = sprintf('FFT Size: %d', N_fft);
summary_text{end+1} = sprintf('Window Length: %.1f s', window_length);
summary_text{end+1} = sprintf('Overlap: %.1f%%', overlap_ratio*100);

text(ax4, 0.05, 0.95, summary_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontSize', 9, 'FontName', 'Courier New', 'BackgroundColor', 'white', 'EdgeColor', 'k');

axis(ax4, 'off');
title(ax4, 'Analysis Summary');

% Save Figure 4
saveas(fig4, fullfile(output_folder, '04_theoretical_analysis.png'));

%% --------------- Save results ---------------
fprintf('Saving frequency analysis results...\n');

% Save analysis data
analysis_results = struct();
analysis_results.meta = struct();
analysis_results.meta.scenario = damage_scenario_selection;
analysis_results.meta.road_class = road_class;
analysis_results.meta.V_kmh = Vehicle_Speed;
analysis_results.meta.fs = Sampling_Freq;
analysis_results.meta.solver_method = solver_method;
analysis_results.meta.timestamp = char(datetime("now","TimeZone","local","Format","yyyy-MM-dd_HH:mm:ss"));

analysis_results.signal_stats = struct();
analysis_results.signal_stats.mean = accel_mean;
analysis_results.signal_stats.std = accel_std;
analysis_results.signal_stats.rms = accel_rms;
analysis_results.signal_stats.max = accel_max;
analysis_results.signal_stats.duration = t_sol(end) - t_sol(1);
analysis_results.signal_stats.samples = length(sprung_mass_accel);

analysis_results.theoretical = struct();
analysis_results.theoretical.body_mode_freq = fn1;
analysis_results.theoretical.wheel_hop_freq = fn2;
analysis_results.theoretical.body_mode_damping = zeta1;
analysis_results.theoretical.wheel_hop_damping = zeta2;

analysis_results.frequency_analysis = struct();
analysis_results.frequency_analysis.dominant_freq = dominant_frequency;
analysis_results.frequency_analysis.peak_frequencies = peak_frequencies;
analysis_results.frequency_analysis.peak_magnitudes = peak_magnitudes;
analysis_results.frequency_analysis.freq_vector = freq_display;
analysis_results.frequency_analysis.fft_magnitude = fft_magnitude_display;
analysis_results.frequency_analysis.fft_phase = fft_phase_display;
analysis_results.frequency_analysis.psd_freq = f_psd_display;
analysis_results.frequency_analysis.psd = psd_display;

analysis_results.analysis_params = struct();
analysis_results.analysis_params.window_length = window_length;
analysis_results.analysis_params.overlap_ratio = overlap_ratio;
analysis_results.analysis_params.max_freq_display = max_freq_display;
analysis_results.analysis_params.min_freq_display = min_freq_display;
analysis_results.analysis_params.fft_size = N_fft;

% Save to file
results_path = fullfile(output_folder, 'frequency_analysis_results.mat');
save(results_path, 'analysis_results', '-v7.3');

% Create summary report
summary_path = fullfile(output_folder, 'frequency_analysis_summary.txt');
fid = fopen(summary_path, 'w');
fprintf(fid, 'FREQUENCY ANALYSIS SUMMARY\n');
fprintf(fid, '=========================\n\n');
fprintf(fid, 'Scenario: %s\n', damage_scenario_selection);
fprintf(fid, 'Road Class: %s\n', road_class);
fprintf(fid, 'Speed: %.1f km/h\n', Vehicle_Speed);
fprintf(fid, 'Sampling Rate: %.0f Hz\n', Sampling_Freq);
fprintf(fid, 'Solver Method: %s\n', solver_method);
fprintf(fid, '\nSIGNAL STATISTICS\n');
fprintf(fid, '=================\n');
fprintf(fid, 'Mean: %.4f m/s²\n', accel_mean);
fprintf(fid, 'Standard Deviation: %.4f m/s²\n', accel_std);
fprintf(fid, 'RMS: %.4f m/s²\n', accel_rms);
fprintf(fid, 'Maximum: %.4f m/s²\n', accel_max);
fprintf(fid, 'Duration: %.2f s\n', t_sol(end) - t_sol(1));
fprintf(fid, 'Samples: %d\n', length(sprung_mass_accel));
fprintf(fid, '\nTHEORETICAL NATURAL FREQUENCIES\n');
fprintf(fid, '===============================\n');
fprintf(fid, 'Body Mode: %.2f Hz (damping: %.3f)\n', fn1, zeta1);
fprintf(fid, 'Wheel Hop: %.2f Hz (damping: %.3f)\n', fn2, zeta2);
fprintf(fid, '\nFREQUENCY ANALYSIS RESULTS\n');
fprintf(fid, '==========================\n');
fprintf(fid, 'Dominant Frequency: %.2f Hz\n', dominant_frequency);
fprintf(fid, 'Body Mode Error: %.2f%%\n', abs(dominant_frequency - fn1)/fn1*100);
fprintf(fid, 'Wheel Hop Error: %.2f%%\n', abs(dominant_frequency - fn2)/fn2*100);
fprintf(fid, '\nDETECTED PEAK FREQUENCIES\n');
fprintf(fid, '=========================\n');
if ~isempty(peak_frequencies)
    for i = 1:length(peak_frequencies)
        fprintf(fid, '%.2f Hz (magnitude: %.3f)\n', peak_frequencies(i), peak_magnitudes(i));
    end
else
    fprintf(fid, 'No significant peaks detected\n');
end
fprintf(fid, '\nANALYSIS PARAMETERS\n');
fprintf(fid, '===================\n');
fprintf(fid, 'Window Length: %.1f s\n', window_length);
fprintf(fid, 'Overlap Ratio: %.1f%%\n', overlap_ratio*100);
fprintf(fid, 'FFT Size: %d\n', N_fft);
fprintf(fid, 'Frequency Range: %.1f - %.1f Hz\n', min_freq_display, max_freq_display);
fprintf(fid, '\nOUTPUT FILES\n');
fprintf(fid, '============\n');
fprintf(fid, 'Results: %s\n', results_path);
fprintf(fid, 'Figure 1: 01_time_domain_analysis.png\n');
fprintf(fid, 'Figure 2: 02_frequency_domain_analysis.png\n');
fprintf(fid, 'Figure 3: 03_spectrogram_analysis.png\n');
fprintf(fid, 'Figure 4: 04_theoretical_analysis.png\n');
fprintf(fid, 'Summary: %s\n', summary_path);
fclose(fid);

fprintf('Frequency analysis completed successfully!\n');
fprintf('Results saved to: %s\n', output_folder);
fprintf('  - Frequency analysis results: %s\n', results_path);
fprintf('  - Figure 1: 01_time_domain_analysis.png\n');
fprintf('  - Figure 2: 02_frequency_domain_analysis.png\n');
fprintf('  - Figure 3: 03_spectrogram_analysis.png\n');
fprintf('  - Figure 4: 04_theoretical_analysis.png\n');
fprintf('  - Summary report: %s\n', summary_path);

%% --------------- Display final statistics ---------------
fprintf('\n========== FREQUENCY ANALYSIS COMPLETED ==========\n');
fprintf('Signal Statistics:\n');
fprintf('  Mean: %.4f m/s²\n', accel_mean);
fprintf('  RMS: %.4f m/s²\n', accel_rms);
fprintf('  Max: %.4f m/s²\n', accel_max);
fprintf('  Duration: %.2f s\n', t_sol(end) - t_sol(1));
fprintf('\nTheoretical Natural Frequencies:\n');
fprintf('  Body Mode: %.2f Hz (damping: %.3f)\n', fn1, zeta1);
fprintf('  Wheel Hop: %.2f Hz (damping: %.3f)\n', fn2, zeta2);
fprintf('\nFrequency Analysis Results:\n');
fprintf('  Dominant Frequency: %.2f Hz\n', dominant_frequency);
fprintf('  Body Mode Error: %.2f%%\n', abs(dominant_frequency - fn1)/fn1*100);
fprintf('  Wheel Hop Error: %.2f%%\n', abs(dominant_frequency - fn2)/fn2*100);
fprintf('  Peak Frequencies Detected: %d\n', length(peak_frequencies));
if ~isempty(peak_frequencies)
    fprintf('  Peak Frequencies:\n');
    for i = 1:min(5, length(peak_frequencies))
        fprintf('    %.2f Hz\n', peak_frequencies(i));
    end
end
fprintf('==================================================\n');