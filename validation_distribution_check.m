
%% VALIDATION: IRI Statistical Distribution (Magnitude Check)
% -------------------------------------------------------------------------
% PURPOSE:
% Verifies that the magnitude of the calculated IRI values follows a
% plausible distribution for real roads. This helps detect calibration
% errors (e.g., wrong suspension parameters causing 2x or 0.5x scaling).
%
% OUTPUTS:
% - Histogram of IRI values with standard road quality bands.
% - Summary statistics (Mean, Median, 95th Percentile).
% -------------------------------------------------------------------------

close all; clear; clc;

%% ========================== LOAD DATA ==========================
iri_file = fullfile('00_Outputs', '05_IRI_Analysis', 'ReconstructedProfile', 'iri_recon_profile.mat');

if ~isfile(iri_file)
    error('IRI file not found: %s', iri_file);
end
data = load(iri_file);
iri_vals = data.iri_rec.segments.IRI_m_per_km;

% Filter out NaNs
iri_vals = iri_vals(~isnan(iri_vals));

%% ========================== STATISTICAL ANALYSIS ==========================
mu    = mean(iri_vals);
med   = median(iri_vals);
p95   = prctile(iri_vals, 95);
sigma = std(iri_vals);

fprintf('\n--- IRI MAGNITUDE REPORT ---\n');
fprintf('Count:  %d segments\n', length(iri_vals));
fprintf('Mean:   %.2f m/km\n', mu);
fprintf('Median: %.2f m/km\n', med);
fprintf('Max:    %.2f m/km\n', max(iri_vals));
fprintf('--------------------------------\n');

% Heuristic Verdicts based on typical paved roads
if med < 0.5
    fprintf('VERDICT: SUSPICIOUSLY LOW. (Are units correct? Is signal lost?)\n');
elseif med > 8.0
    fprintf('VERDICT: SUSPICIOUSLY HIGH. (Check suspension parameters or units)\n');
elseif med > 1.5 && med < 4.0
    fprintf('VERDICT: PLAUSIBLE (Typical Highway/City Road range)\n');
else
    fprintf('VERDICT: ATYPICAL (Might be unpaved or damaged road)\n');
end

%% ========================== VISUALIZATION ==========================
figure('Name', 'IRI Distribution Check', 'Color', 'w', 'Position', [150 150 800 500]);
hold on; grid on; box on;

% 1. Draw Standard Road Class Bands (Background)
% Green: Good (<2.5), Yellow: Fair (2.5-4.0), Red: Poor (>4.0)
max_y = length(iri_vals) * 0.4; % Estimate histogram height for shading scaling
yl = [0, max_y * 1.5]; % Y-limits

patch([0 2.5 2.5 0], [yl(1) yl(1) yl(2) yl(2)], [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Good (<2.5)');
patch([2.5 4.0 4.0 2.5], [yl(1) yl(1) yl(2) yl(2)], [1 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Fair (2.5-4.0)');
patch([4.0 20 20 4.0], [yl(1) yl(1) yl(2) yl(2)], [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Poor (>4.0)');

% 2. Histogram
histogram(iri_vals, 30, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.8, 'EdgeColor', 'w', 'DisplayName', 'Calculated IRI');

% 3. Mark Median
xline(med, 'k--', ['Median: ' num2str(med, '%.2f')], 'LineWidth', 2, 'DisplayName', 'Median');

xlim([0, max(10, p95*1.2)]); % Dynamic X limit based on data, but at least 10
ylim(yl);
xlabel('IRI (m/km)');
ylabel('Count (Number of Segments)');
title('IRI Magnitude Distribution vs. Standard Classes');
legend('Location', 'northeast');

fprintf('Plot generated.\n');