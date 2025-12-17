
%% VALIDATION: PSD Slope Verification (Spectral Consistency)
% -------------------------------------------------------------------------
% PURPOSE:
% Validates the physical realism of a reconstructed road profile by 
% analyzing its frequency content. According to ISO 8608, a natural road
% profile must have a PSD slope of approximately -2 (log-log scale).
%
%   - Slope approx -2.0 : Valid Road Profile (Physical)
%   - Slope approx  0.0 : White Noise (Sensor/Solver Noise)
%   - Slope approx < -3.0 : Integration Drift (Low-frequency instability)
%
% INPUTS:
%   - 'recon_spatial.mat' (contains 'sp' struct with x and zr)
%
% OUTPUTS:
%   - Diagnostic plot (PSD vs ISO Reference).
%   - Console report with the calculated slope and verdict.
% -------------------------------------------------------------------------

close all; clear; clc;

%% 1) Load Reconstructed Profile
% -------------------------------------------------------------------------
% Path to the reconstructed spatial profile
% Check common locations for the file
file_candidates = { ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'SpatialMapping', 'recon_spatial.mat'), ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'PathB_FrequencyDomain', 'recon_spatial.mat'), ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'PathA_TimeDomain', 'recon_spatial.mat'), ...
    'recon_spatial.mat'
};

mat_path = '';
for k = 1:numel(file_candidates)
    if exist(file_candidates{k}, 'file')
        mat_path = file_candidates{k};
        break;
    end
end

if isempty(mat_path)
    error('Could not find "recon_spatial.mat". Run the reconstruction first.');
end

fprintf('Loading profile from: %s\n', mat_path);
data = load(mat_path);

% Extract Data
if isfield(data, 'sp')
    sp = data.sp;
elseif isfield(data, 'pkg') && isfield(data.pkg, 'spatial')
    sp = data.pkg.spatial; % Legacy support
else
    error('Invalid .mat structure. Expected struct "sp" or "pkg.spatial".');
end

% Choose profile to analyze (Filtered is usually preferred for final check)
if isfield(sp, 'zr_x_filt')
    zr = sp.zr_x_filt;
    name_str = 'Filtered Profile';
elseif isfield(sp, 'zr_x_raw')
    zr = sp.zr_x_raw;
    name_str = 'Raw Profile';
else
    % Legacy field names
    zr = sp.zr_filt;
    name_str = 'Filtered Profile';
end

% Extract Spatial Sampling Frequency
if isfield(sp, 'fs_spatial')
    fs_spatial = sp.fs_spatial; % samples per meter
else
    dx = mean(diff(sp.x));
    fs_spatial = 1 / dx;
end

% Ensure column vector and remove DC offset (mean)
zr = zr(:);
zr = detrend(zr, 'constant'); 

%% 2) Calculate Power Spectral Density (PSD)
% -------------------------------------------------------------------------
% Use Welch's method for a clean, smooth spectrum
nfft = 2^nextpow2(length(zr));
window = hamming(floor(length(zr)/8)); % 8 segments window
noverlap = floor(length(window)/2);

[pxx, n_spatial] = pwelch(zr, window, noverlap, nfft, fs_spatial);

% pxx units: m^3 (displacement PSD)
% n_spatial units: cycles/m (spatial frequency)

%% 3) Analyze Slope (Linear Fit in Log-Log Domain)
% -------------------------------------------------------------------------
% We only fit the "Road Band" of interest, typically 0.05 to 2.0 cycles/m.
% This avoids the very low freq (drift) and very high freq (sensor noise floor).

fit_band = [0.05, 2.0]; % cycles/m (Wavelengths 20m down to 0.5m)

% Find indices for the fit
idx_fit = find(n_spatial >= fit_band(1) & n_spatial <= fit_band(2));

if isempty(idx_fit)
    warning('Spatial sampling rate too low to analyze road roughness band. Using all available data.');
    idx_fit = 2:length(n_spatial); % Fallback to all (skip DC)
end

% Prepare Log Data for fitting
log_n = log10(n_spatial(idx_fit));
log_p = log10(pxx(idx_fit));

% Linear Regression: log(P) = m * log(n) + c
p_poly = polyfit(log_n, log_p, 1);
slope = p_poly(1);
intercept = p_poly(2);

% Generate best-fit line for plotting
fit_y = polyval(p_poly, log_n);
fit_p = 10.^fit_y;
fit_n = n_spatial(idx_fit);

%% 4) Interpretation & Verdict
% -------------------------------------------------------------------------
fprintf('\n--- SPECTRAL CONSISTENCY REPORT ---\n');
fprintf('Calculated Slope (w): %.2f\n', slope);

verdict = '';
color_verdict = '';

if slope > -0.5
    verdict = 'FAIL: WHITE NOISE (Slope ~ 0)';
    detail = 'The profile is dominated by sensor noise. IRI will be falsely high.';
    color_verdict = 'r';
elseif slope < -3.5
    verdict = 'FAIL: DRIFT / INTEGRATION ERROR (Slope < -3)';
    detail = 'Low-frequency drift is dominating. Profile is wandering vertically.';
    color_verdict = 'm';
elseif slope >= -3.0 && slope <= -1.5
    verdict = 'PASS: VALID ROAD PROFILE (Slope ~ -2)';
    detail = 'The spectral shape matches physical road topology (ISO 8608).';
    color_verdict = 'g';
else
    verdict = 'WARNING: ATYPICAL PROFILE';
    detail = 'Slope is somewhat consistent but outside typical ISO bounds.';
    color_verdict = [0.9 0.6 0]; % Orange
end

fprintf('Verdict: %s\n', verdict);
fprintf('Details: %s\n', detail);

%% 5) Visualization
% -------------------------------------------------------------------------
figure('Name', 'Spectral Consistency Check', 'Color', 'w', 'Position', [100 100 900 600]);

% Main Log-Log Plot
loglog(n_spatial, pxx, 'k', 'LineWidth', 1.5, 'DisplayName', 'Reconstructed PSD');
hold on;
grid on;

% Plot the Fit Line
loglog(fit_n, fit_p, 'r--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Fit Slope = %.2f', slope));

% Plot Reference Slope (-2) for visual comparison
% Anchor the reference line to the middle of the fit data to make it visible nearby
ref_mid_n = median(fit_n);
ref_mid_p = median(fit_p);
ref_y = ref_mid_p * (fit_n ./ ref_mid_n).^(-2);
loglog(fit_n, ref_y, 'b:', 'LineWidth', 1.5, 'DisplayName', 'Ideal Road (Slope -2)');

% Formatting
xlabel('Spatial Frequency n (cycles/m)');
ylabel('Displacement PSD G_d(n) (m^3)');
title(['Spectral Consistency Check: ' name_str]);
legend('Location', 'best');
xlim([0.01, max(n_spatial)]);

% Add Text Box with Verdict
dim = [0.15 0.35 0.3 0.1];
str = {['Slope w = ' num2str(slope, '%.2f')], verdict};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', 'EdgeColor', color_verdict, 'LineWidth', 2, 'FontSize', 11);

% Highlight the Fit Band
xline(fit_band(1), 'g--', 'HandleVisibility', 'off');
xline(fit_band(2), 'g--', 'HandleVisibility', 'off');
text(fit_band(1), max(pxx)/100, ' Analysis Band Start', 'Rotation', 90, 'VerticalAlignment', 'bottom');

fprintf('\nPlot generated.\n');