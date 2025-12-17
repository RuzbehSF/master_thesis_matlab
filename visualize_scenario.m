%% 02) Scenario & 3D Visualization (Post-processing)
% -------------------------------------------------------------------------
% PURPOSE:
%   Visualize the selected damage scenario after running the main
%   quarter-car simulation script.
%
%   Figure 1:
%     - Subplot 1: Generated road profile (elevation vs distance)
%     - Subplot 2: Sprung vertical acceleration vs time
%     - Subplot 3: Damage Scenario Layout (schematic, along x)
%
%   Figure 2:
%     - 3D view of the road surface + overlaid damage scenario cues
%
% REQUIREMENTS:
%   Run the main simulation script first so that:
%       00_Outputs/01_Simulation/SimulationData/simulation_results.mat
%   exists and contains 'results'.
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1) Load simulation results

results_path = fullfile('00_Outputs','01_Simulation','SimulationData','simulation_results.mat');

if ~isfile(results_path)
    error('Could not find results file:\n  %s\nRun the main simulation script first.', results_path);
end

S = load(results_path);
if ~isfield(S, 'results')
    error('File "%s" does not contain a struct named "results".', results_path);
end

results = S.results;

%% 2) Extract core data

% Spatial road profile
x  = results.spatial_road(:);     % m
zr = results.zr.x(:);             % road profile with damages (m)

if isfield(results.zr, 'base_x') && ~isempty(results.zr.base_x)
    zr_base = results.zr.base_x(:);   % ISO baseline (no deterministic damages)
else
    zr_base = nan(size(zr));
end

if isempty(x)
    error('results.spatial_road is empty. Nothing to visualize.');
end

L = x(end) - x(1);                % approximate road length (m)

% Time & acceleration
t = results.time(:);              % s
if isfield(results.outputs, 'accel_clean')
    az_clean = results.outputs.accel_clean(:);   % m/s^2 (sprung vertical accel, clean)
else
    error('results.outputs.accel_clean not found – cannot plot sprung vertical acceleration.');
end

if isfield(results.outputs, 'accel_meas')
    az_meas = results.outputs.accel_meas(:);     % noisy channel, if present
else
    az_meas = [];
end

% Meta
scenario   = string(results.meta.scenario);    % e.g., "Bumps_Only"
road_class = string(results.meta.road_class);  % e.g., "A"
V_kmh      = results.meta.V_kmh;
fs         = results.meta.fs;

%% 3) Styling

set(0,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',11,...
       'DefaultAxesLineWidth',0.8,'DefaultLineLineWidth',1.25,...
       'DefaultFigureColor','w','DefaultTextInterpreter','none');

COL.road     = [0.20 0.45 0.85];
COL.base     = [0.55 0.55 0.55];
COL.acc_clean= [0.00 0.00 0.00];
COL.acc_meas = [0.85 0.33 0.10];

COL.bandLow  = [0.85 0.92 1.00];
COL.bandMed  = [0.90 0.95 1.00];
COL.bandHi   = [0.80 0.88 0.98];
COL.bandExt  = [0.75 0.82 0.96];

COL.bump     = [0.00 0.45 0.74];
COL.pothole  = [0.85 0.33 0.10];
COL.rut      = [0.47 0.67 0.19];
COL.crack    = [0.55 0.00 0.55];
COL.corr     = [0.93 0.69 0.13];
COL.edge     = [0.30 0.30 0.30];

symLim = @(y) max(1e-6, max(abs(y)))*1.1;

%% 4) Figure 1: 3 Subplots (road, acceleration, layout)

fig1 = figure('Name','Scenario Overview','Position',[80 60 1200 900]);
tl1 = tiledlayout(fig1,3,1,'TileSpacing','compact','Padding','compact');

title(tl1, sprintf('Scenario: %s   |   ISO Class: %s   |   V=%.0f km/h, L=%.0f m, fs=%.0f Hz', ...
    strrep(scenario,'_',' '), road_class, V_kmh, L, fs), ...
    'FontSize',12,'FontWeight','bold');

%% 4.1 Subplot 1: Generated road profile

ax1 = nexttile(tl1,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on'); grid(ax1,'minor');

y_mm = zr*1000;
plot(ax1, x, y_mm, 'Color', COL.road, 'DisplayName','Road (with damages)');

if all(isfinite(zr_base)) && any(abs(zr_base) > 1e-12)
    plot(ax1, x, zr_base*1000, '--', 'Color', COL.base, 'DisplayName','ISO baseline');
end

xlabel(ax1,'Distance x (m)');
ylabel(ax1,'Elevation (mm)');
ylim(ax1, symLim(y_mm)*[-1 1]);
legend(ax1,'Location','best','Box','off');
title(ax1,'Generated Road Profile');

%% 4.2 Subplot 2: Sprung vertical acceleration

ax2 = nexttile(tl1,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on'); grid(ax2,'minor');

plot(ax2, t, az_clean, 'Color', COL.acc_clean, 'DisplayName','Sprung accel (clean)');

if ~isempty(az_meas) && numel(az_meas) == numel(t)
    plot(ax2, t, az_meas, '-', 'Color', COL.acc_meas, 'DisplayName','Sprung accel (measured)');
end

xlabel(ax2,'Time t (s)');
ylabel(ax2,'a_z (m/s^2)');
title(ax2,'Sprung Vertical Acceleration');
legend(ax2,'Location','best','Box','off');

%% 4.3 Subplot 3: Damage Scenario Layout

ax3 = nexttile(tl1,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
xlabel(ax3,'Distance x (m)');
set(ax3,'YTick',[]);
ylim(ax3,[0 1]);
xlim(ax3,[x(1) x(end)]);

drawScenarioLayout(ax3, scenario, x, L, COL);
title(ax3,'Damage Scenario Layout (schematic)');

fprintf('\nFigure 1 created: road, sprung acceleration, and layout.\n');

%% 5) Figure 2: 3D visualization of road + damage scenario

fig2 = figure('Name','3D Road & Damage Scenario','Position',[150 100 1200 800]);
ax3d = axes(fig2); hold(ax3d,'on'); box(ax3d,'on'); grid(ax3d,'on');

% Build a simple lane-width extrusion of the 1D profile
laneWidth = 3.5;            % m, arbitrary lane width for visualization
y_lane = linspace(-laneWidth/2, laneWidth/2, 8);

[Xs, Ys] = meshgrid(x, y_lane);
Zs = repmat(zr.', numel(y_lane), 1);   % same profile across lane width

surf(ax3d, Xs, Ys, Zs*1000);  % convert to mm for Z
shading(ax3d,'interp');
colormap(ax3d, 'parula');
alpha(ax3d,0.95);

xlabel(ax3d,'Distance x (m)');
ylabel(ax3d,'Lateral y (m)');
zlabel(ax3d,'Elevation (mm)');
title(ax3d, sprintf('3D Road Surface | Scenario: %s', strrep(scenario,'_',' ')));

view(ax3d, [-35 25]);
axis(ax3d,'tight');
grid(ax3d,'on');

% Overlay damage scenario cues in 3D
drawScenario3D(ax3d, scenario, x, zr, laneWidth, COL);

camlight(ax3d,'headlight'); lighting(ax3d,'gouraud');

fprintf('Figure 2 created: 3D road & scenario visualization.\n\n');

%% ========================================================================
%% Local helper functions
%% ========================================================================

function drawScenarioLayout(ax, scenario, x, L, COL)
% Draw a 2D schematic layout of the selected scenario on axes 'ax'.

    s = lower(strtrim(scenario));
    s = strrep(s,' ','_');

    sevNames = {'Low','Medium','High','Extreme'};

    hold(ax,'on');

    switch s
        case 'smooth'
            text(mean(x),0.5,'Smooth (ISO roughness only, no deterministic damages)', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');

        case 'highway_smooth'
            drawHighwaySmooth_2D(ax, L, COL);

        case 'city_urban'
            drawCityUrban_2D(ax, L, COL);

        case 'progressive_roughness'
            drawProgressiveRoughness_2D(ax, x, COL);

        % ===== NEW SCENARIOS (2D LAYOUT) =================================
        case 'body_resonance_track'
            % Sinusoidal track tuned to body mode (schematic band)
            width = 0.6*L;
            center = 0.5*L;
            drawBand2D(ax, center, width, 0.35, 0.65, COL.bandMed);
            text(center, 0.5, 'Body resonance sinusoidal track', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');

        case 'wheelhop_resonance_track'
            % High frequency sinusoidal track (wheel-hop)
            width = 0.6*L;
            center = 0.5*L;
            drawBand2D(ax, center, width, 0.35, 0.65, COL.bandHi);
            text(center, 0.5, 'Wheel-hop resonance track', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');

        case 'cobblestone_street'
            % Cobblestone-like section in middle of road
            xs = 0.20*L;
            xe = 0.70*L;
            drawBand2D(ax, (xs+xe)/2, xe-xs, 0.30, 0.70, COL.bandMed);
            text((xs+xe)/2, 0.5, 'Cobblestone section', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');

        case 'gravel_road'
            % Long rough gravel segment
            drawBand2D(ax, 0.5*L, L, 0.25, 0.75, COL.bandExt);
            text(0.5*L, 0.5, 'Gravel road (rough, loose surface)', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');

        case 'rail_crossing'
            % Localized crossing in middle
            xc = 0.5*L;
            drawMarker2D(ax, xc, 0.5, 'x', COL.corr, 'Rail / tram crossing');
            xs = 0.5*L - 0.15*L;
            xe = 0.5*L + 0.15*L;
            drawBand2D(ax, (xs+xe)/2, xe-xs, 0.25, 0.35, COL.bandLow);

        case 'random_damage_field'
            % Random mix across most of the length (schematic band)
            xs = 0.10*L;
            xe = 0.90*L;
            drawBand2D(ax, (xs+xe)/2, xe-xs, 0.20, 0.80, COL.bandMed);
            text((xs+xe)/2, 0.5, 'Random mix of bumps, potholes, ruts, etc.', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');

        % ===== ORIGINAL SCENARIOS ========================================
        case 'bumps_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            for i=1:4
                drawMarker2D(ax, centers(i), 0.25, '^', COL.bump, ...
                    sprintf('Bump %s Severity', sevNames{i}));
            end

        case 'potholes_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            for i=1:4
                drawMarker2D(ax, centers(i), 0.5, 'v', COL.pothole, ...
                    sprintf('Pothole %s Severity', sevNames{i}));
            end

        case 'rutting_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            yRow = 0.6;
            for i=1:4
                drawBand2D(ax, centers(i), L/10, yRow-0.08, yRow+0.08, COL.bandMed);
                text(centers(i), yRow, sprintf('Rut %s Severity', sevNames{i}), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
            end

        case 'crack_joints_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            for i=1:4
                drawMarker2D(ax, centers(i), 0.7, '|', COL.crack, ...
                    sprintf('Crack/Joint %s Severity', sevNames{i}));
            end

        case 'block_cracking_only'
            bands = bands4(L);
            drawSeverityBands2D(ax, bands, sevNames, COL);
            textMidBands2D(ax, bands, {'Block cracking (Low)', 'Block cracking (Medium)', ...
                                       'Block cracking (High)', 'Block cracking (Extreme)'});

        case 'fatigue_crocodile_only'
            bands = bands4(L);
            drawSeverityBands2D(ax, bands, sevNames, COL);
            textMidBands2D(ax, bands, {'Fatigue / crocodile (Low)', 'Fatigue / crocodile (Medium)', ...
                                       'Fatigue / crocodile (High)', 'Fatigue / crocodile (Extreme)'});

        case 'reflection_cracking_only'
            bands = bands4(L);
            drawSeverityBands2D(ax, bands, sevNames, COL);
            textMidBands2D(ax, bands, {'Reflection cracking (Low)', 'Reflection cracking (Medium)', ...
                                       'Reflection cracking (High)', 'Reflection cracking (Extreme)'});

        case 'corrugation_shoving_only'
            bands = bands4(L);
            drawSeverityBands2D(ax, bands, sevNames, COL);
            textMidBands2D(ax, bands, {'Corrugation + shoving (Low)', 'Corrugation + shoving (Medium)', ...
                                       'Corrugation + shoving (High)', 'Corrugation + shoving (Extreme)'});

        case 'depression_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            for i=1:4
                drawMarker2D(ax, centers(i), 0.5, 'v', COL.rut, ...
                    sprintf('Depression %s Severity', sevNames{i}));
            end

        case 'swelling_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            for i=1:4
                drawMarker2D(ax, centers(i), 0.5, '^', COL.rut, ...
                    sprintf('Swelling %s Severity', sevNames{i}));
            end

        case 'edge_fatigue_only'
            bands = bands4(L);
            drawSeverityBands2D(ax, bands, sevNames, COL);
            textMidBands2D(ax, bands, {'Edge fatigue (Low)', 'Edge fatigue (Medium)', ...
                                       'Edge fatigue (High)', 'Edge fatigue (Extreme)'});

        case {'low_severity','medium_severity','high_severity','extreme_severity'}
            sev = scenarioSeverity(s);
            drawBand2D(ax, 0.5*L, L, 0.2, 0.8, severityColor(sev, COL));
            text(0.5*L, 0.5, sprintf('All damage types (%s severity)', sev), ...
                'HorizontalAlignment','center','VerticalAlignment','middle');

        case 'all_damages_mixed_severity'
            bands = bands4(L);
            sevOrder = {'Low','Medium','High','Extreme'};
            for i=1:4
                drawBand2D(ax, mean(bands(i,:)), diff(bands(i,:)), 0.2, 0.8, severityColor(sevOrder{i}, COL));
                text(mean(bands(i,:)), 0.5, sprintf('All damages (%s)', sevOrder{i}), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
            end

        otherwise
            text(mean(x),0.5, sprintf('No dedicated layout for "%s".\nRoad profile above is still valid.', scenario), ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

function drawScenario3D(ax, scenario, x, zr, laneWidth, COL)
% Overlay 3D cues on top of the extruded road surface.

    s = lower(strtrim(scenario));
    s = strrep(s,' ','_');

    hold(ax,'on');

    switch s
        case 'smooth'
            return;

        case 'highway_smooth'
            drawHighwaySmooth_3D(ax, x, zr, laneWidth, COL);

        case 'city_urban'
            drawCityUrban_3D(ax, x, zr, laneWidth, COL);

        case 'progressive_roughness'
            drawProgressiveRoughness_3D(ax, x, zr, laneWidth, COL);

        % ===== NEW SCENARIOS (3D) ========================================
        case 'body_resonance_track'
            L = x(end) - x(1);
            xs = 0.20*L + x(1);
            xe = 0.80*L + x(1);
            draw3DBand(ax, xs, xe, laneWidth, zr, x, COL.bandMed);

        case 'wheelhop_resonance_track'
            L = x(end) - x(1);
            xs = 0.20*L + x(1);
            xe = 0.80*L + x(1);
            draw3DBand(ax, xs, xe, laneWidth, zr, x, COL.bandHi);

        case 'cobblestone_street'
            L = x(end) - x(1);
            xs = 0.20*L + x(1);
            xe = 0.70*L + x(1);
            draw3DBand(ax, xs, xe, laneWidth, zr, x, COL.bandMed);

        case 'gravel_road'
            L = x(end) - x(1);
            xs = 0.05*L + x(1);
            xe = 0.95*L + x(1);
            draw3DBand(ax, xs, xe, laneWidth, zr, x, COL.bandExt);

        case 'rail_crossing'
            L = x(end) - x(1);
            xc = 0.5*L + x(1);
            xs = xc - 0.15*L;
            xe = xc + 0.15*L;
            draw3DBand(ax, xs, xe, laneWidth, zr, x, COL.bandLow);
            draw3DMarkerLine(ax, xc, laneWidth, zr, x, COL.corr);

        case 'random_damage_field'
            L = x(end) - x(1);
            xs = 0.10*L + x(1);
            xe = 0.90*L + x(1);
            draw3DBand(ax, xs, xe, laneWidth, zr, x, COL.bandMed);

        % ===== ORIGINAL SCENARIOS ========================================
        case 'bumps_only'
            centers = linspace(0.15*(x(end)-x(1)), 0.85*(x(end)-x(1)), 4);
            centers = centers + x(1);
            for i=1:numel(centers)
                draw3DMarkerLine(ax, centers(i), laneWidth, zr, x, COL.bump);
            end

        case 'potholes_only'
            centers = linspace(0.15*(x(end)-x(1)), 0.85*(x(end)-x(1)), 4);
            centers = centers + x(1);
            for i=1:numel(centers)
                draw3DMarkerLine(ax, centers(i), laneWidth, zr, x, COL.pothole);
            end

        case 'rutting_only'
            bands = bands4(x(end)-x(1));
            bands = bands + x(1);
            for i=1:4
                draw3DBand(ax, bands(i,1), bands(i,2), laneWidth, zr, x, COL.rut);
            end

        case 'crack_joints_only'
            centers = linspace(0.15*(x(end)-x(1)), 0.85*(x(end)-x(1)), 4);
            centers = centers + x(1);
            for i=1:numel(centers)
                draw3DMarkerLine(ax, centers(i), laneWidth, zr, x, COL.crack);
            end

        case {'block_cracking_only','fatigue_crocodile_only',...
              'reflection_cracking_only','corrugation_shoving_only',...
              'edge_fatigue_only'}
            bands = bands4(x(end)-x(1));
            bands = bands + x(1);
            sevColors = {COL.bandLow, COL.bandMed, COL.bandHi, COL.bandExt};
            for i=1:4
                draw3DBand(ax, bands(i,1), bands(i,2), laneWidth, zr, x, sevColors{i});
            end

        case {'depression_only','swelling_only'}
            centers = linspace(0.15*(x(end)-x(1)), 0.85*(x(end)-x(1)), 4);
            centers = centers + x(1);
            ccol = COL.rut;
            for i=1:numel(centers)
                draw3DMarkerLine(ax, centers(i), laneWidth, zr, x, ccol);
            end

        case {'low_severity','medium_severity','high_severity','extreme_severity'}
            sev = scenarioSeverity(s);
            c = severityColor(sev, COL);
            draw3DBand(ax, x(1), x(end), laneWidth, zr, x, c);

        case 'all_damages_mixed_severity'
            bands = bands4(x(end)-x(1));
            bands = bands + x(1);
            sevColors = {COL.bandLow, COL.bandMed, COL.bandHi, COL.bandExt};
            for i=1:4
                draw3DBand(ax, bands(i,1), bands(i,2), laneWidth, zr, x, sevColors{i});
            end

        otherwise
            return;
    end
end

%% -------- geometry helpers (2D) ---------

function bands = bands4(L)
% Four equal-length bands along [0, L].
    bands = [0       L/4; ...
             L/4     L/2; ...
             L/2   3*L/4; ...
             3*L/4   L];
end

function drawBand2D(ax, center, width, y0, y1, colorVal)
    xs = center - width/2;
    xe = center + width/2;
    patch(ax, [xs xe xe xs], [y0 y0 y1 y1], colorVal, ...
        'EdgeColor','none','FaceAlpha',0.4);
end

function drawMarker2D(ax, x0, y0, mark, colorVal, labelStr)
    % Vertical guide line
    line(ax, [x0 x0], [0 y0], 'Color', colorVal, 'LineStyle',':');

    % Arrow marker
    plot(ax, x0, y0, mark, 'MarkerSize',8, ...
        'MarkerFaceColor', colorVal, 'MarkerEdgeColor', colorVal);

    % Text label above the arrow (e.g., "Pothole Low Severity")
    text(ax, x0, y0 + 0.06, labelStr, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');
end

function drawSeverityBands2D(ax, bands, sevNames, COL)
    for i=1:4
        xs = bands(i,1);
        xe = bands(i,2);
        c = severityColor(sevNames{i}, COL);
        drawBand2D(ax, (xs+xe)/2, xe-xs, 0.25, 0.75, c);
    end
end

function textMidBands2D(ax, bands, labels)
    for i=1:4
        xs = bands(i,1);
        xe = bands(i,2);
        text(ax, (xs+xe)/2, 0.5, labels{i}, ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

function drawHighwaySmooth_2D(ax, L, COL)
    joints = 0:50:L;
    joints = joints(joints>0 & joints<L);
    utility = [0.2 0.4 0.6 0.8] * L;

    for j = joints
        drawMarker2D(ax, j, 0.4, '|', COL.crack, 'Expansion joint');
    end
    for u = utility
        drawMarker2D(ax, u, 0.7, 's', COL.pothole, 'Utility cut / patch');
    end
end

function drawCityUrban_2D(ax, L, COL)
    speed_bumps = [0.15 0.35 0.55 0.75 0.95] * L;
    potholes    = [0.10 0.25 0.45 0.65 0.85] * L;
    util_cuts   = [0.20 0.40 0.60 0.80] * L;
    old_sec     = [0.30 0.70] * L;

    for sb = speed_bumps
        drawMarker2D(ax, sb, 0.2, '^', COL.bump, 'Speed bump');
    end
    for p = potholes
        drawMarker2D(ax, p, 0.5, 'v', COL.pothole, 'Pothole');
    end
    for u = util_cuts
        drawMarker2D(ax, u, 0.75, 's', COL.rut, 'Utility cut / patch');
    end
    for s = old_sec
        len = 0.1*L;
        xs = s;
        xe = min(s+len,L);
        drawBand2D(ax, (xs+xe)/2, xe-xs, 0.65, 0.90, COL.bandHi);
        text((xs+xe)/2, 0.8, 'Block cracking zone', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

function drawProgressiveRoughness_2D(ax, x, COL)
    L = x(end) - x(1);
    edges = linspace(x(1), x(end), 9);
    labelsAH = 'ABCDEFGH';

    for k = 1:8
        xs = edges(k);
        xe = edges(k+1);
        drawBand2D(ax, (xs+xe)/2, xe-xs, 0.2, 0.8, COL.bandLow);
        text((xs+xe)/2, 0.5, labelsAH(k), ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
    end
end

%% -------- geometry helpers (3D) ---------

function draw3DMarkerLine(ax, x0, laneWidth, zr, x, colorVal)
    % Draw a vertical "cut" line across the lane at given x0
    z0 = interp1(x, zr, x0, 'linear', 'extrap');
    yv = [-laneWidth/2, laneWidth/2];
    zv = [z0 z0]*1000 + 2;  % little offset above surface (mm)
    plot3(ax, [x0 x0], yv, zv, 'LineWidth',2, 'Color', colorVal);
end

function draw3DBand(ax, xs, xe, laneWidth, zr, x, colorVal)
    % Draw a translucent overlay band on the surface between xs and xe
    idx = x>=xs & x<=xe;
    if ~any(idx)
        return;
    end
    xseg = [xs xe xe xs];
    yseg = [-laneWidth/2 -laneWidth/2 laneWidth/2 laneWidth/2];
    zbase = interp1(x, zr, [xs xe xe xs], 'linear', 'extrap')*1000;
    zseg = zbase + 1.5;   % offset slightly above road (mm)
    patch(ax, xseg, yseg, zseg, colorVal, ...
        'FaceAlpha',0.35,'EdgeColor','none');
end

function drawHighwaySmooth_3D(ax, x, zr, laneWidth, COL)
    joints = 0:50:(x(end)-x(1));
    joints = joints(joints>0 & joints<(x(end)-x(1)));
    joints = joints + x(1);
    utility = [0.2 0.4 0.6 0.8]*(x(end)-x(1));
    utility = utility + x(1);

    for j = joints
        draw3DMarkerLine(ax, j, laneWidth, zr, x, COL.crack);
    end
    for u = utility
        draw3DMarkerLine(ax, u, laneWidth, zr, x, COL.pothole);
    end
end

function drawCityUrban_3D(ax, x, zr, laneWidth, COL)
    L = x(end) - x(1);

    speed_bumps = [0.15 0.35 0.55 0.75 0.95]*L + x(1);
    potholes    = [0.10 0.25 0.45 0.65 0.85]*L + x(1);
    util_cuts   = [0.20 0.40 0.60 0.80]*L + x(1);
    old_sec     = [0.30 0.70]*L + x(1);    % just for bands

    for sb = speed_bumps
        draw3DMarkerLine(ax, sb, laneWidth, zr, x, COL.bump);
    end
    for p = potholes
        draw3DMarkerLine(ax, p, laneWidth, zr, x, COL.pothole);
    end
    for u = util_cuts
        draw3DMarkerLine(ax, u, laneWidth, zr, x, COL.rut);
    end
    for s = old_sec
        len = 0.1*L;
        xs = s;
        xe = min(s+len,x(end));
        draw3DBand(ax, xs, xe, laneWidth, zr, x, COL.bandHi);
    end
end

function drawProgressiveRoughness_3D(ax, x, zr, laneWidth, COL)
    L = x(end) - x(1);
    edges = linspace(x(1), x(end), 9);
    for k = 1:8
        xs = edges(k);
        xe = edges(k+1);
        draw3DBand(ax, xs, xe, laneWidth, zr, x, COL.bandLow);
    end
end

%% -------- misc helpers ---------

function c = severityColor(sev, COL)
    sev = lower(sev);
    switch sev
        case 'low'
            c = COL.bandLow;
        case 'medium'
            c = COL.bandMed;
        case 'high'
            c = COL.bandHi;
        case 'extreme'
            c = COL.bandExt;
        otherwise
            c = [0.9 0.9 0.9];
    end
end

function sev = scenarioSeverity(s)
    if contains(s,'low','IgnoreCase',true)
        sev = 'Low';
    elseif contains(s,'medium','IgnoreCase',true)
        sev = 'Medium';
    elseif contains(s,'high','IgnoreCase',true)
        sev = 'High';
    elseif contains(s,'extreme','IgnoreCase',true)
        sev = 'Extreme';
    else
        sev = 'Unknown';
    end
end
