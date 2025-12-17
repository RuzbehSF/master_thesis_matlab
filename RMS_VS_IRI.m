%% VALIDATION: RMS "Sanity Check" (Linear Correlation)
% -------------------------------------------------------------------------
% PURPOSE:
% Validates the physical consistency of the reconstructed IRI by correlating
% it with the Root Mean Square (RMS) of the raw vertical acceleration signal.
%
% LOGIC:
% 1. Loads the raw acceleration signal (input to the reconstruction).
% 2. Loads the pre-calculated per-segment IRI values.
% 3. Calculates the RMS of the acceleration for those EXACT same segments.
% 4. Plots IRI vs. RMS Accel and calculates the R-squared correlation.
%
% EXPECTED RESULT:
% - A strong linear correlation (R^2 > 0.8).
% - If R^2 is low (< 0.5), the reconstruction is likely failing or noisy.
% - If linear but offset, a calibration factor is needed.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ========================== HYPERPARAMETERS ==========================
% (1) Paths to Input Files
sim_file = fullfile('00_Outputs', '01_Simulation', 'SimulationData', 'simulation_results.mat');
iri_file = fullfile('00_Outputs', '05_IRI_Analysis', 'ReconstructedProfile', 'iri_recon_profile.mat');

% (2) Settings
accel_source        = 'meas';  % 'meas' (measured/noisy) or 'clean' (ground truth/clean)
min_segment_points  = 20;      % Minimum number of accel points required to calculate RMS for a segment

% (3) Output
save_plot           = true;
plot_path           = fullfile('00_Outputs', '04_Validation', 'RMS_Check_Plot.png');

%% ========================== 1. LOAD DATA ==========================
fprintf('Loading data...\n');

% Load Simulation Data (Acceleration & Distance)
if ~isfile(sim_file)
    error('Simulation results file not found: %s', sim_file);
end
sim_data = load(sim_file);
results  = sim_data.results;

% Extract Acceleration
if strcmp(accel_source, 'meas')
    accel = results.outputs.accel_meas(:);
    fprintf('Using MEASURED (Noisy) Acceleration.\n');
else
    accel = results.outputs.accel_clean(:);
    fprintf('Using CLEAN (Ground Truth) Acceleration.\n');
end

% Extract Distance Vector (Ensure it matches accel length)
if isfield(results, 'x_m') && length(results.x_m) == length(accel)
    x_dist = results.x_m(:);
elseif isfield(results, 'spatial_road')
    x_dist = results.spatial_road(:);
    if length(x_dist) ~= length(accel)
        V = results.meta.V_kmh / 3.6;
        t = results.time(:);
        x_dist = t * V;
    end
else
    error('Could not determine distance vector matching acceleration data.');
end

% Load IRI Data
if ~isfile(iri_file)
    error('IRI results file not found: %s', iri_file);
end
iri_data = load(iri_file);
iri_rec  = iri_data.iri_rec;

% Extract Segment Table
if isfield(iri_rec, 'segments') && istable(iri_rec.segments)
    SegData    = iri_rec.segments;
    segs_start = SegData.x_start_m;
    segs_end   = SegData.x_end_m;
    segs_iri   = SegData.iri_m_per_km;
elseif isfield(iri_rec.segments, 'edges_m')
    edges      = iri_rec.segments.edges_m; % Nx2 matrix [start, end]
    segs_start = edges(:,1);
    segs_end   = edges(:,2);
    segs_iri   = iri_rec.segments.IRI_m_per_km;
else
    error('Unknown IRI segments structure.');
end

%% ========================== 2. CALCULATE RMS PER SEGMENT ==========================
fprintf('Calculating RMS Acceleration for %d segments...\n', length(segs_iri));

num_segs = length(segs_iri);
segs_rms = nan(num_segs, 1);

for i = 1:num_segs
    x_s = segs_start(i);
    x_e = segs_end(i);

    idx = find(x_dist >= x_s & x_dist < x_e);

    if length(idx) >= min_segment_points
        segment_accel = accel(idx);
        segs_rms(i)   = rms(segment_accel);    % Standard RMS (includes DC)
    else
        segs_rms(i)   = NaN;
    end
end

valid_mask = ~isnan(segs_rms) & ~isnan(segs_iri);
Y_IRI      = segs_iri(valid_mask);
X_RMS      = segs_rms(valid_mask);

if isempty(Y_IRI)
    error('No valid segments found overlapping between Distance and Acceleration data.');
end

%% ========================== 3. STATISTICAL VALIDATION ==========================
% Linear Regression: IRI = a * RMS + b
mdl     = fitlm(X_RMS, Y_IRI);
R2      = mdl.Rsquared.Ordinary;
Coeffs  = mdl.Coefficients.Estimate; % [Intercept; Slope]
Slope   = Coeffs(2);
Intercept = Coeffs(1);

fprintf('\n--- VALIDATION RESULTS ---\n');
fprintf('Number of Segments: %d\n', length(Y_IRI));
fprintf('Correlation (R^2) : %.4f\n', R2);
fprintf('Linear Fit        : IRI = %.2f * RMS + %.2f\n', Slope, Intercept);

if R2 > 0.8
    fprintf('STATUS: PASS (Strong correlation)\n');
elseif R2 > 0.5
    fprintf('STATUS: WARNING (Moderate correlation)\n');
else
    fprintf('STATUS: FAIL (Weak correlation - check reconstruction or alignment)\n');
end

%% ========================== 4. PLOTTING (ENHANCED) ==========================
% ---- Global style ----
baseFont = 11;
lineLW   = 1.6;

set(0,'defaultAxesFontName','Helvetica', ...
      'defaultAxesFontSize',baseFont, ...
      'defaultLineLineWidth',lineLW, ...
      'defaultAxesBox','off', ...
      'defaultAxesXGrid','on', ...
      'defaultAxesYGrid','on', ...
      'defaultAxesXMinorGrid','on', ...
      'defaultAxesYMinorGrid','on');

% Color palette (soft, colorblind-friendly-ish)
C.blue  = [0.00 0.45 0.74];
C.red   = [0.80 0.25 0.25];
C.gray  = [0.55 0.55 0.60];
C.green = [0.20 0.62 0.20];

fig = figure('Name','RMS Sanity Check', ...
             'Color','w', ...
             'Units','pixels', ...
             'Position',[100 100 900 640]);
fig.Renderer = 'painters';

tlo = tiledlayout(fig,1,1,'TileSpacing','compact','Padding','compact');
ax  = nexttile(tlo); hold(ax,'on');

% Scatter: use slight jitter and alpha to avoid overplotting
markerSize = 45;
scatter(ax, X_RMS, Y_IRI, markerSize, ...
    'MarkerFaceColor',C.blue, ...
    'MarkerEdgeColor',[0 0 0], ...
    'MarkerFaceAlpha',0.85, ...
    'MarkerEdgeAlpha',0.6, ...
    'DisplayName','Segments');

% Fit line and 95% confidence band
x_fit = linspace(min(X_RMS), max(X_RMS), 200).';
[y_pred, y_ci] = predict(mdl, x_fit);

% Confidence band (shaded area)
fill(ax, [x_fit; flipud(x_fit)], ...
         [y_ci(:,1); flipud(y_ci(:,2))], ...
         C.blue, ...
         'FaceAlpha',0.10, ...
         'EdgeColor','none', ...
         'DisplayName','95% CI');

% Regression line
plot(ax, x_fit, y_pred, '-', ...
    'Color',C.red, ...
    'LineWidth',2, ...
    'DisplayName',sprintf('Linear fit (R^2 = %.2f)', R2));

% Optional: dashed horizontal line at mean IRI (visual cue)
y_mean = mean(Y_IRI,'omitnan');
plot(ax, [min(X_RMS) max(X_RMS)], [y_mean y_mean], ':', ...
    'Color',C.gray, ...
    'LineWidth',1.0, ...
    'HandleVisibility','off');

% Axes labels and title
xlabel(ax,'RMS Vertical Acceleration (m/s^2)','FontSize',baseFont+1);
ylabel(ax,'Calculated IRI (m/km)','FontSize',baseFont+1);
title(ax,'RMS Sanity Check: IRI vs Vertical Acceleration','FontSize',baseFont+3);

% Make axes a bit tighter with some margin
axis(ax,'tight');
xlim(ax, [min(X_RMS)*0.98, max(X_RMS)*1.02]);
ylim(ax, [min(Y_IRI)*0.98, max(Y_IRI)*1.05]);

% Nicely styled legend
leg = legend(ax,'Location','northwest');
set(leg,'Box','off','FontSize',baseFont-1,'ItemTokenSize',[12 8]);

% Text box with fit equation, placed inside the axes
txt_str = sprintf('IRI = %.2f · RMS + %.2f\nR^2 = %.3f', Slope, Intercept, R2);
x_txt   = 0.97*min(xlim(ax)) + 0.03*max(xlim(ax));
y_txt   = 0.97*max(ylim(ax)) - 0.13*range(ylim(ax));
text(ax, x_txt, y_txt, txt_str, ...
    'FontSize',baseFont, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'BackgroundColor',[1 1 1 0.85], ...
    'Margin',5, ...
    'EdgeColor',[0.8 0.8 0.8]);

% Make sure grids look good
set(ax,'LineWidth',0.9, ...
       'GridAlpha',0.25, ...
       'MinorGridAlpha',0.12);

% Save Plot (high-quality PNG)
if save_plot
    save_dir = fileparts(plot_path);
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end

    [~,~,ext] = fileparts(plot_path);
    if isempty(ext)
        plot_path = [plot_path '.png'];
    end

    try
        exportgraphics(fig, plot_path, 'Resolution', 300);
    catch
        % Fallback for older MATLAB
        print(fig, plot_path, '-dpng', '-r300');
    end
    fprintf('Plot saved to: %s\n', plot_path);
end
