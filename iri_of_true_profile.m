
%% 10) Calculate IRI on the Original (Ground Truth) Profile
% -------------------------------------------------------------------------
% PURPOSE:
% This script calculates the International Roughness Index (IRI) for the
% original, ground-truth road profile. This serves as the benchmark value
% for the road's roughness.
%
% CORE LOGIC (How it works):
% 1.  Loads the configuration file ('recon_cfg.mat') to access the ground
%     truth road profile (x_true, zr_true).
%
% 2.  Simulates the standard "Golden Car" quarter-car model driving over
%     the true road profile at the IRI standard speed of 80 km/h. The
%     Golden Car has a precisely defined set of parameters used for all
%     official IRI calculations.
%
% 3.  The simulation calculates the vehicle's dynamic response, specifically
%     the relative velocity between the sprung mass (body) and unsprung
%     mass (wheel).
%
% 4.  It calculates the overall IRI by taking the integrated absolute value
%     of this relative velocity, normalized by the length of the road,
%     according to the standard IRI definition.
%
% 5.  It also calculates IRI values for fixed-length segments (e.g., every
%     100 meters) to show how roughness varies along the road.
%
% 6.  The IRI values are classified into categories (e.g., "Very Good",
%     "Good", "Fair") and saved to .mat and .csv files.
%
% INPUTS:
% - '02_Reconstruction Inputs/recon_cfg.mat': For the ground truth profile.
%
% OUTPUTS:
% - '10_IRI True Profile/iri_original_profile.mat': Contains the calculated
%   IRI values and vehicle response.
% - '10_IRI True Profile/iri_original_segments.csv': A table of IRI values
%   per road segment.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ============================== HYPERPARAMETERS ==============================
% (1) Input discovery (new location first, legacy fallback)
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};

% (2) This script’s output folder
out_root = fullfile('00_Outputs', '05_IRI_Analysis', 'TrueProfile');
fig_dir  = fullfile(out_root,'figs');

% (3) Golden Car parameters (canonical defaults)
GC.ms_corner = 250;  % kg    (sprung mass, per wheel)
GC.mu = 35;          % kg    (unsprung mass)
GC.ks = 63.3e3;      % N/m   (suspension stiffness)
GC.cs = 6.53e3;      % Ns/m  (suspension damping)
GC.kt = 653e3;       % N/m   (tire stiffness)
GC.ct = 0;           % Ns/m  (tire damping; IRI standard ~0)

% (4) IRI simulation speed & numerics
iri_speed_kmh     = 80;     % km/h (IRI is defined at 80 km/h)
target_fs_time    = [];     % Hz; [] -> use native sampling from profile at Viri
warmup_length_m   = 0;      % meters ignored at the beginning for metrics

% (5) Segmentation for reporting
segment_length_m  = 10;     % meters per segment for IRI reporting (typical 10/20/50/100 m)

% (6) Units to report and plot
%   'm_per_km' -> report in m/km only
%   'mm_per_m' -> report in mm/m only
%   'both'     -> save both in MAT/CSV; plots use m/km
iri_units_report  = 'm_per_km';

% (7) Plots & export
make_plots = true;
export_png = true;

%% ------------------------------- Load inputs --------------------------------
% Pick first available cfg
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file')==2, cfg_path = cfg_candidates{k}; break; end
end
assert(~isempty(cfg_path), 'Missing recon_cfg.mat in expected locations.');

C = load(cfg_path); cfg = C.cfg;

if ~exist(out_root,'dir')
    mkdir(out_root);
end
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Ground truth profile (required)
assert(isfield(cfg,'ground') && isfield(cfg.ground,'x') && ~isempty(cfg.ground.x), ...
    'Ground-truth profile not found in cfg.ground.x. Please run your simulator to create it.');
x_true  = cfg.ground.x(:);        % distance [m]
zr_true = cfg.ground.zr(:);       % elevation [m]

% Ensure strictly increasing x (monotonic for interp/ODE timing)
[x_true, uniq_idx] = unique(x_true, 'stable');
zr_true = zr_true(uniq_idx);

% IRI speed & timing
Viri   = iri_speed_kmh / 3.6;            % m/s
t_true = (x_true - x_true(1)) / Viri;    % time mapping at IRI speed

% Optional resampling in time (usually not necessary if dx is uniform)
if ~isempty(target_fs_time)
    t_end = t_true(end);
    t     = (0:1/target_fs_time:t_end).';
    zr    = interp1(t_true, zr_true, t, 'linear', 'extrap');
else
    t  = t_true;
    zr = zr_true;
end

% Derivative of road elevation w.r.t. time
% 250 mm Moving Average Filter (Simulate Tire Enveloping)
% Calculate effective spatial sampling interval (dx = V * dt)
dt_avg = mean(diff(t));
dx_val = Viri * dt_avg;          % meters per sample
window_sz = round(0.250 / dx_val); % 250mm window size

if window_sz > 1
    % Apply standard moving average (boxcar filter)
    zr = movmean(zr, window_sz);
    fprintf('Applied 250mm moving average filter (Window: %d samples, dx=%.4fm)\n', window_sz, dx_val);
end

% Road time derivative (calculated on the now-smoothed profile)
dzr_dt = deriv1(zr, t);

%% --------------------------- Golden Car simulation --------------------------
[A,B] = quarter_AB(GC);
Cout  = eye(4); Dout = zeros(4, size(B,2));
sys   = ss(A,B,Cout,Dout);    % continuous LTI

% Inputs U = [zr, dzr_dt]
U = [zr, dzr_dt];

% Simulate with zero initial conditions
X = lsim(sys, U, t, zeros(4,1));

% States mapping (consistent with quarter_AB): [ys, dys, yus, dyus]
ys   = X(:,1);
dys  = X(:,2);
yus  = X(:,3);
dyus = X(:,4);

% Relative suspension velocity
v_rel = dys - dyus;    % z_ṡ - z_u̇

% Distance array aligned with t (IRI speed)
x = t * Viri + x_true(1);

% Warmup mask
m_warm = x >= (x(1) + warmup_length_m);

% Total path length (excluding warmup) for normalization
L_total = max(x(m_warm)) - min(x(m_warm));
L_total = max(L_total, eps);

% Overall IRI [m/km] (base unit for classification)
IRI_total_m_per_km = (1000 / L_total) * trapz(t(m_warm), abs(v_rel(m_warm)));

%% -------------------- Segment-wise IRI on fixed-length bins -----------------
segs = make_segments(x, segment_length_m);
nSeg = size(segs,1);
IRI_seg_m_per_km = nan(nSeg,1);
RCI_seg          = nan(nSeg,1);

for k = 1:nSeg
    xs = segs(k,1); xe = segs(k,2);
    m  = (x >= xs) & (x < xe);
    if any(m)
        Lk = xe - xs; Lk = max(Lk, eps);
        IRI_seg_m_per_km(k) = (1000 / Lk) * trapz(t(m), abs(v_rel(m)));
    end
end

% Segment centers for plotting
seg_centers = mean(segs,2);

%% ------------------------- Units & Qualification ----------------------------
% Convert to requested units (numerically, 1 m/km == 1 mm/m)
[overall_val, overall_units_lbl] = iri_convert_units(IRI_total_m_per_km, iri_units_report);
seg_vals = iri_convert_units(IRI_seg_m_per_km, iri_units_report);

% Classify overall and segments using base m/km thresholds
[overall_class, overall_color] = iri_classify_single(IRI_total_m_per_km);
[seg_class, seg_color, seg_class_idx] = iri_classify(IRI_seg_m_per_km);

% Also keep both-unit versions in outputs for convenience
overall_mm_per_m = IRI_total_m_per_km;   % numerically identical
seg_mm_per_m     = IRI_seg_m_per_km;

%% ------------------------------- Save outputs -------------------------------
iri = struct();
iri.info.model       = 'Golden Car (quarter-car)';
iri.info.params      = GC;
iri.info.speed_kmh   = iri_speed_kmh;
iri.info.warmup_m    = warmup_length_m;
iri.info.seg_len_m   = segment_length_m;
iri.info.units_mode  = iri_units_report;

iri.series.x_m       = x;
iri.series.zr_m      = zr;
iri.series.v_rel_mps = v_rel;
iri.series.t_s       = t;

iri.overall.IRI_m_per_km  = IRI_total_m_per_km;
iri.overall.IRI_mm_per_m  = overall_mm_per_m;
iri.overall.value_display = overall_val;
iri.overall.units_display = overall_units_lbl;
iri.overall.class_name    = overall_class;
iri.overall.color_rgb     = overall_color;

iri.segments.edges_m        = segs;
iri.segments.centers_m      = seg_centers;
iri.segments.IRI_m_per_km   = IRI_seg_m_per_km;
iri.segments.IRI_mm_per_m   = seg_mm_per_m;
iri.segments.value_display  = seg_vals;
iri.segments.class_name     = seg_class;
iri.segments.class_index    = seg_class_idx;
iri.segments.color_rgb      = seg_color;

save(fullfile(out_root,'iri_original_profile.mat'), 'iri','-v7.3');
fprintf('Saved: %s  (overall IRI = %.3f m/km, class = %s)\n', ...
    fullfile(out_root,'iri_original_profile.mat'), IRI_total_m_per_km, overall_class);

% CSV for segments (include both units and class)
Tseg = table((1:nSeg).', segs(:,1), segs(:,2), (segs(:,2)-segs(:,1)), ...
             IRI_seg_m_per_km, seg_mm_per_m, seg_vals(:), seg_class(:), seg_class_idx(:), ...
    'VariableNames', {'segment_idx','x_start_m','x_end_m','length_m', ...
                      'iri_m_per_km','iri_mm_per_m','iri_value_display','iri_class','class_idx'});
writetable(Tseg, fullfile(out_root,'iri_original_segments.csv'));
fprintf('Saved: %s\n', fullfile(out_root,'iri_original_segments.csv'));

%% --------------------------------- Plots -----------------------------------
if make_plots
    % Fig 1 — Original profile with segment grids
    f1 = mkfig('IRI (original) — Profile & segments', [100 70 1100 540]);
    plot(x, zr, 'k', 'LineWidth', 1.2); grid on; box on; formatAxes();
    title('Original road profile  z_r(x)');
    xlabel('Distance (m)'); ylabel('Elevation (m)');
    yl = ylim; hold on;
    for k = 1:nSeg
        plot([segs(k,1) segs(k,1)], yl, ':', 'Color',[0.75 0.75 0.75]);
    end
    if export_png
        exportgraphics(f1, fullfile(fig_dir,'23_iri_original_profile.png'), 'Resolution', 180);
    end

    % Fig 2 — Segment IRI bars with background bands (no legend, no y-limit forcing)
    f2 = mkfig('IRI (original) — Segment bars (qualified)', [130 90 1100 560]);
    ax = gca; hold(ax,'on'); box(ax,'on'); formatAxes(ax);

    % Bars colored by class (in display units)
    b = bar(ax, seg_centers, seg_vals, 1.0, 'EdgeColor','none');
    b.FaceColor = 'flat';
    b.CData     = seg_color;   % Nx3 per-segment class colors

    % Overall line (already in display units)
    yline(ax, overall_val, 'k--', ...
        sprintf('overall = %.3f %s (%s)', overall_val, overall_units_lbl, overall_class), ...
        'LabelHorizontalAlignment','left');

    % Draw bands AFTER plotting (uses current auto y-limits; does not change them)
    draw_iri_bands_auto(ax, overall_units_lbl, 0.12);

    title(ax, sprintf('Segment IRI (%d m) — qualified background', segment_length_m));
    xlabel(ax, 'Segment center (m)'); ylabel(ax, sprintf('IRI (%s)', overall_units_lbl));
    if export_png
        exportgraphics(f2, fullfile(fig_dir,'24_iri_original_segments_banded.png'), 'Resolution', 180);
    end

    % Fig 3 — Running IRI vs distance (plots in display units)
    f3 = mkfig('IRI (original) — Running IRI', [160 110 980 520]);
    running_mpkm = running_iri_per_km(t, x, v_rel, warmup_length_m);
    running_disp = iri_convert_units(running_mpkm, iri_units_report);
    plot(x, running_disp, 'm', 'LineWidth', 1.2); grid on; box on; formatAxes();
    title('Running IRI (from start)'); xlabel('Distance (m)'); ylabel(sprintf('IRI (%s)', overall_units_lbl));
    yline(overall_val, 'k--', 'final overall');
    if export_png
        exportgraphics(f3, fullfile(fig_dir,'25_iri_original_running.png'), 'Resolution', 180);
    end

    % Fig 4 — Relative suspension velocity trace (warmup highlighted)
    f4 = mkfig('IRI (original) — Relative suspension velocity', [190 130 1100 520]);
    plot(x, v_rel, 'Color',[0.25 0.45 0.85], 'LineWidth', 1.0); grid on; box on; formatAxes(); hold on;
    title('Relative suspension velocity  v_{rel}(x) = z_ṡ - z_u̇');
    xlabel('Distance (m)'); ylabel('Velocity (m/s)');
    xline(x(1)+warmup_length_m, 'r--', sprintf('warmup = %dm', warmup_length_m));
    if export_png
        exportgraphics(f4, fullfile(fig_dir,'26_iri_original_vrel.png'), 'Resolution', 180);
    end
end


%% =============================== FUNCTIONS ===============================
function [A,B] = quarter_AB(p)
% Quarter-car (Golden Car) continuous-time state-space:
% state x = [ys; dys; yus; dyus], inputs u = [zr; dzr]
%   ddys  = [ -cs(dys-dyus) - ks(ys-yus) ] / ms
%   ddyus = [  cs(dys-dyus) + ks(ys-yus) - ct(dyus-dzr) - kt(yus-zr) ] / mu
    ms = p.ms_corner;
    mu = p.mu;
    ks = p.ks;
    cs = p.cs;
    kt = p.kt;
    ct = p.ct;

    A = [0    1     0       0;
         -ks/ms -cs/ms  ks/ms   cs/ms;
          0    0     0       1;
          ks/mu cs/mu -(ks+kt)/mu  -(cs+ct)/mu];
    B = [0  0;
         0  0;
         0  0;
         kt/mu  ct/mu];
end

function dz = deriv1(z,t)
% First derivative using safe forward/central differences on nonuniform t.
    z = z(:); t = t(:);
    n = numel(z);
    dz = zeros(n,1);
    if n < 2
        return;
    end
    dt = diff(t);
    dt(dt==0) = eps;
    dz(1)      = (z(2)-z(1)) / dt(1);
    dz(2:n-1)  = (z(3:n) - z(1:n-2)) ./ (t(3:n) - t(1:n-2));
    dz(n)      = (z(n)-z(n-1)) / dt(end);
end

function segs = make_segments(x, Lseg)
% Build [start, end] edges covering x-range with fixed length Lseg (last may be shorter).
    xmin = x(1); xmax = x(end);
    edges = xmin : Lseg : (xmax + 1e-9);
    if edges(end) < xmax
        edges = [edges xmax];
    end
    segs = [edges(1:end-1).' edges(2:end).'];
end

function iri_running = running_iri_per_km(t, x, v_rel, warmup_m)
% Compute cumulative IRI along the path; ignore first warmup_m meters.
    iri_running = nan(size(x));
    x0 = x(1) + warmup_m;
    m0 = x >= x0;
    if ~any(m0)
        iri_running(:) = 0;
        return;
    end
    idx0 = find(m0,1,'first');
    cumI = cumtrapz(t(idx0:end), abs(v_rel(idx0:end)));
    Ls   = x(idx0:end) - x(idx0);
    Ls(Ls<=0) = eps;
    iri_running(idx0:end) = 1000 * (cumI ./ Ls);  % m/km
    iri_running(1:idx0-1) = NaN;
end

function [val_disp, units_lbl] = iri_convert_units(val_m_per_km, mode)
% Convert IRI values to requested display units.
% Note: numerically, 1 m/km == 1 mm/m.
    switch lower(mode)
        case 'm_per_km'
            val_disp = val_m_per_km;
            units_lbl = 'm/km';
        case 'mm_per_m'
            val_disp = val_m_per_km;   % same numeric
            units_lbl = 'mm/m';
        case 'both'
            val_disp = val_m_per_km;   % plots will use m/km
            units_lbl = 'm/km';
        otherwise
            error('Unknown iri_units_report mode: %s', mode);
    end
end

function [name, rgb] = iri_classify_single(iri_mpkm)
% Classify a single IRI (in m/km) and return class name + color.
    [names, colors, edges] = iri_classes_def();
    if iri_mpkm < edges(2)
        idx = 1;
    elseif iri_mpkm < edges(3)
        idx = 2;
    elseif iri_mpkm < edges(4)
        idx = 3;
    elseif iri_mpkm <= edges(5)
        idx = 4;
    else
        idx = 5;
    end
    name = names{idx};
    rgb  = colors(idx,:);
end

function [names, colors, idx] = iri_classify(iri_mpkm_vec)
% Vectorized classifier for segment arrays.
    [namesLUT, colorsLUT, edges] = iri_classes_def();
    v = iri_mpkm_vec(:);
    idx = zeros(size(v));
    for i=1:numel(v)
        if v(i) < edges(2)
            idx(i) = 1;
        elseif v(i) < edges(3)
            idx(i) = 2;
        elseif v(i) < edges(4)
            idx(i) = 3;
        elseif v(i) <= edges(5)
            idx(i) = 4;
        else
            idx(i) = 5;
        end
    end
    names  = cellfun(@(k) namesLUT{k}, num2cell(idx), 'UniformOutput', false);
    colors = colorsLUT(idx,:);   % Nx3 per-segment colors
end

function [names, colors, edges] = iri_classes_def()
% Table: Very Good <2, Good 2–4, Fair 4–6, Poor 6–10, Very Poor >10  (m/km)
    names  = {'Very Good','Good','Fair','Poor','Very Poor'};
    colors = [  0.00 0.60 0.05;   % green
                0.70 0.82 0.20;   % yellow-green
                0.95 0.70 0.20;   % orange
                0.90 0.00 0.00;   % red
                0.60 0.10 0.10];  % dark red
    edges  = [0 2 4 6 10];        % helper thresholds (m/km)
end

function hBands = draw_iri_bands_auto(ax, units_lbl, alpha_bg)
% Draw semi-transparent IRI class bands WITHOUT changing y-limits.
% Call AFTER plotting bars/lines so the axis has auto-scaled.
% Returns Patch handles (already placed behind other content).
    if nargin < 3 || isempty(alpha_bg), alpha_bg = 0.12; end
    [~, baseColors, edges] = iri_classes_def();        % edges = [0 2 4 6 10]
    colors = 1 - (1 - baseColors) * 0.35;              % lighten for backgrounds

    % Use current limits (do NOT modify them)
    xl = xlim(ax); yl = ylim(ax);

    % Upper for each band: [2 4 6 10  (current top)]
    uppers = [edges(2:end) yl(2)];
    lowers = edges;

    hBands = gobjects(5,1);
    for i = 1:5
        y0 = max(yl(1), lowers(i));
        y1 = min(yl(2), uppers(i));
        if y1 > y0
            hBands(i) = patch(ax, ...
                [xl(1) xl(2) xl(2) xl(1)], [y0 y0 y1 y1], colors(i,:), ...
                'EdgeColor','none', 'FaceAlpha',alpha_bg, ...
                'HandleVisibility','off', 'HitTest','off');
        end
    end

    % Threshold lines only if visible
    for e = edges(2:end)
        if e > yl(1) && e < yl(2)
            yline(ax, e, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75, 'HandleVisibility','off');
        end
    end

    % Compact labels at the right edge
    DX = diff(xl);
    x_text = xl(2) - 0.01*DX;
    labels = {sprintf('< 2.0 %s', units_lbl), sprintf('2.0–4.0 %s', units_lbl), ...
              sprintf('4.0–6.0 %s', units_lbl), sprintf('6.0–10.0 %s', units_lbl), ...
              sprintf('> 10.0 %s', units_lbl)};
    for i = 1:5
        % Only place label if the band exists
        y0 = max(yl(1), lowers(i));
        y1 = min(yl(2), uppers(i));
        if y1 > y0
            y_mid = (y0 + y1)/2;
            text(ax, x_text, y_mid, labels{i}, ...
                 'HorizontalAlignment','right','VerticalAlignment','middle', ...
                 'FontSize',8, 'Color',[0.25 0.25 0.25], 'HitTest','off');
        end
    end

    % Keep bands behind everything else (only valid handles)
    hb = hBands(isgraphics(hBands));
    if ~isempty(hb)
        uistack(hb, 'bottom');
    end
    ax.Layer = 'top';
end

function f = mkfig(name, pos)
    f = figure('Name', name, 'Color', 'w', 'Position', pos);
end

function formatAxes(ax)
    if nargin<1, ax = gca; end
    ax.LineWidth = 1.0;
    ax.FontName  = 'Calibri';
    ax.FontSize  = 11;
    grid(ax,'on'); box(ax,'on');
end