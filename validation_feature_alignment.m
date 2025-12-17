%% VALIDATION: Feature Detection (Event Alignment)
% -------------------------------------------------------------------------
% PURPOSE:
% Verifies that high-roughness events in the calculated IRI profile align
% spatially with high-energy events in the raw acceleration signal.
%
% This detects:
% 1. Spatial Drift: If peaks consistently drift apart, your speed/distance
%    calculation is wrong.
% 2. Filtering Lag: If IRI peaks consistently lag behind acceleration,
%    your smoothing filters have too much phase delay.
%
% OUTPUTS:
% - Normalized overlay plot of Acceleration Envelope vs. IRI.
% - Peak matching report (distance error between events).
% -------------------------------------------------------------------------

close all; clear; clc;

%% ========================== HYPERPARAMETERS ==========================
% (1) Input Files
sim_file = fullfile('00_Outputs', '01_Simulation', 'SimulationData', 'simulation_results.mat');
iri_file = fullfile('00_Outputs', '05_IRI_Analysis', 'ReconstructedProfile', 'iri_recon_profile.mat');

% (2) Settings
window_size_m    = 10;    % Window size for smoothing acceleration (match IRI segment size)
peak_prominence  = 0.2;   % Min prominence (0-1 scale) to consider a "major event"
match_tolerance_m = 20;   % Max distance allowed between accel peak and IRI peak to be a "match"

%% ========================== 1. LOAD & PREP DATA ==========================
fprintf('Loading data...\n');
sim = load(sim_file);
iri = load(iri_file);

% Extract Raw Acceleration
accel = sim.results.outputs.accel_meas(:);

% Extract Distance for Acceleration (Time * Speed)
t_sim = sim.results.time(:);
V_mps = sim.results.meta.V_kmh / 3.6;
x_accel = t_sim * V_mps;

% Extract IRI Profile (segment centers)
iri_vals    = iri.iri_rec.segments.IRI_m_per_km;
iri_centers = iri.iri_rec.segments.centers_m;

%% ========================== 2. PROCESS SIGNALS ==========================
% 1. Acceleration "Energy Envelope"
fs = sim.results.meta.fs;
samples_per_window = max(1, round(window_size_m / (V_mps/fs)));
accel_envelope = movmean(abs(accel), samples_per_window);

% 2. Interpolate acceleration envelope onto IRI centers
accel_interp = interp1(x_accel, accel_envelope, iri_centers, 'linear', 'extrap');

% 3. Normalize both signals to [0,1]
norm_accel = (accel_interp - min(accel_interp)) / range(accel_interp);
norm_iri   = (iri_vals   - min(iri_vals))   / range(iri_vals);

%% ========================== 3. PEAK MATCHING ==========================
[pks_A, locs_A] = findpeaks(norm_accel, iri_centers, 'MinPeakProminence', peak_prominence);
[pks_I, locs_I] = findpeaks(norm_iri,   iri_centers, 'MinPeakProminence', peak_prominence);

fprintf('\n--- EVENT ALIGNMENT REPORT ---\n');
fprintf('Found %d Acceleration Events and %d IRI Events.\n', length(locs_A), length(locs_I));

matched_diffs = [];
match_pairs   = []; % [x_accel_peak, x_iri_peak]

for i = 1:length(locs_A)
    [min_dist, idx] = min(abs(locs_I - locs_A(i)));
    if min_dist <= match_tolerance_m
        diff_val = locs_I(idx) - locs_A(i);  % >0 => IRI lags behind
        matched_diffs = [matched_diffs; diff_val]; %#ok<AGROW>
        match_pairs   = [match_pairs; locs_A(i), locs_I(idx)]; %#ok<AGROW>
        fprintf('Event @ %.1fm: IRI offset = %+.2fm\n', locs_A(i), diff_val);
    end
end

if isempty(matched_diffs)
    avg_shift = NaN;
    status_str = 'FAIL (No matching events found)';
    fprintf('WARNING: No matching events found! Check synchronization.\n');
else
    avg_shift = mean(matched_diffs);
    fprintf('------------------------------------------------\n');
    fprintf('AVERAGE SPATIAL SHIFT: %+.2f meters\n', avg_shift);
    if abs(avg_shift) > 5
        status_str = 'FAIL (Significant drift/lag detected)';
        fprintf('STATUS: FAIL (Significant drift/lag detected)\n');
    else
        status_str = 'PASS (Good spatial alignment)';
        fprintf('STATUS: PASS (Good spatial alignment)\n');
    end
end

%% ========================== 4. VISUALIZATION (ENHANCED) ==========================
% Global style
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

% Color palette
C.accel  = [0.00 0.45 0.74];   % blue
C.iri    = [0.80 0.25 0.25];   % red
C.gray   = [0.55 0.55 0.60];
C.green  = [0.20 0.62 0.20];
C.orange = [0.90 0.60 0.00];

fig = figure('Name','Feature Alignment Check', ...
             'Color','w', ...
             'Units','pixels', ...
             'Position',[100 100 1000 600]);
fig.Renderer = 'painters';

tlo = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

%% 4a) Top panel: signals, peaks, and matches
ax1 = nexttile(tlo,1); hold(ax1,'on');

% Slightly smoothed accel envelope for readability (data unchanged)
win_smooth = max(3, round(numel(norm_accel)/400));
norm_accel_smooth = smooth(norm_accel, win_smooth);

% Main traces
h_accel = plot(ax1, iri_centers, norm_accel_smooth, '-', ...
    'Color',C.accel, ...
    'LineWidth',1.4, ...
    'DisplayName','Accel energy (input)');

h_iri = plot(ax1, iri_centers, norm_iri, '-', ...
    'Color',C.iri, ...
    'LineWidth',1.6, ...
    'DisplayName','Calculated IRI (output)');

% Shaded bands around matched events (visual alignment aid)
yl = [0 1.02];    % tight y-range for normalized signals
set(ax1,'YLim',yl);
for k = 1:size(match_pairs,1)
    xc = mean(match_pairs(k,:));          % center between accel & IRI peak
    w  = match_tolerance_m;               % half-width of shaded band
    patch(ax1, ...
        [xc-w xc+w xc+w xc-w], ...
        [yl(1) yl(1) yl(2) yl(2)], ...
        [0.92 0.92 0.92], ...
        'FaceAlpha',0.25, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');
end

% Bring lines above patches
uistack(h_accel,'top');
uistack(h_iri,'top');

% Peaks (different markers for clarity)
if ~isempty(locs_A)
    plot(ax1, locs_A, pks_A, 'o', ...
        'MarkerFaceColor',C.accel, ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',6, ...
        'DisplayName','Accel peaks');
end
if ~isempty(locs_I)
    plot(ax1, locs_I, pks_I, 's', ...
        'MarkerFaceColor',C.iri, ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',6, ...
        'DisplayName','IRI peaks');
end

% Match connection lines (muted so they don’t dominate)
for k = 1:size(match_pairs,1)
    xa = match_pairs(k,1);
    xi = match_pairs(k,2);

    % Get corresponding amplitudes (nearest neighbors if exact match missing)
    [~,ia] = min(abs(iri_centers - xa));
    [~,ii] = min(abs(iri_centers - xi));
    ya = norm_accel(ia);
    yi = norm_iri(ii);

    plot(ax1, [xa xi], [ya yi], '-', ...
        'Color',[0.20 0.62 0.20], ...
        'LineWidth',1.3, ...
        'HandleVisibility','off');
end

% Axes cosmetics
xlim(ax1, [min(iri_centers) max(iri_centers)]);
ylim(ax1, [-0.05 1.05]);

xlabel(ax1,'Distance (m)','FontSize',baseFont+1);
ylabel(ax1,'Normalized magnitude (0–1)','FontSize',baseFont+1);
title(ax1, sprintf('Feature alignment (avg shift = %+.2f m)', avg_shift), ...
      'FontSize',baseFont+3,'FontWeight','bold');

leg1 = legend(ax1,'Location','northwest');
set(leg1,'Box','off','FontSize',baseFont-1,'ItemTokenSize',[12 8]);

set(ax1,'LineWidth',0.9, ...
        'GridAlpha',0.25, ...
        'MinorGridAlpha',0.12);

% Summary box anchored in figure coordinates (never overlaps curves)
summary_str = sprintf( ...
    'Accel events: %d\nIRI events:   %d\nMatches:      %d\nAvg shift:    %+.2f m\n%s', ...
    length(locs_A), length(locs_I), length(matched_diffs), avg_shift, status_str);

annotation(fig,'textbox', ...
    'Units','normalized', ...
    'Position',[0.15 0.78 0.18 0.18], ...  % top-left corner
    'String',summary_str, ...
    'FontSize',baseFont, ...
    'BackgroundColor',[1 1 1 0.85], ...
    'Margin',6, ...
    'EdgeColor',[0.85 0.85 0.85], ...
    'LineWidth',1.0);

%% 4b) Bottom panel: histogram of spatial offsets
ax2 = nexttile(tlo,2); hold(ax2,'on');

if ~isempty(matched_diffs)
    % Symmetric binning around zero
    span = max(abs(matched_diffs));
    span = max(span, match_tolerance_m);        % at least tolerance span
    bin_w = max(0.5, span/6);                   % reasonable bin width
    edges = (-span-bin_w):bin_w:(span+bin_w);

    histogram(ax2, matched_diffs, edges, ...
        'FaceColor',C.accel, ...
        'FaceAlpha',0.75, ...
        'EdgeColor','none');

    y_l = ylim(ax2);

    % Zero reference
    xline(ax2, 0, '-', ...
        'Color',C.gray, ...
        'LineWidth',1.1, ...
        'HandleVisibility','off');

    % Average shift
    xline(ax2, avg_shift, '--', ...
        'Color',C.iri, ...
        'LineWidth',2.0, ...
        'DisplayName','Average shift');

else
    text(ax2,0.5,0.5,'No matched events to plot.', ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color',C.gray, ...
        'FontSize',baseFont+1);
end

xlabel(ax2,'Spatial offset (IRI peak – accel peak) [m]','FontSize',baseFont+1);
ylabel(ax2,'Count','FontSize',baseFont+1);
title(ax2,'Distribution of spatial offsets between matched events', ...
      'FontSize',baseFont+2);

set(ax2,'LineWidth',0.9, ...
        'GridAlpha',0.25, ...
        'MinorGridAlpha',0.12);

fprintf('\nPlot generated.\n');
