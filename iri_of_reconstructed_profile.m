
%% 11) Calculate IRI on the Reconstructed Profile
% -------------------------------------------------------------------------
% PURPOSE:
% This script calculates the International Roughness Index (IRI) for the
% final, reconstructed road profile. The result can be compared to the IRI
% of the true profile to evaluate the performance of the reconstruction in
% terms of this industry-standard metric.
%
% CORE LOGIC (How it works):
% 1.  Loads the final reconstructed spatial profile from
%     'recon_spatial.mat'. The user can choose to use either the "raw" or
%     the "filtered" version of the profile.
%
% 2.  Just like the previous script, it simulates the standard "Golden Car"
%     quarter-car model driving over this reconstructed road profile at the
%     IRI standard speed of 80 km/h.
%
% 3.  The simulation calculates the relative velocity between the Golden
%     Car's sprung and unsprung masses.
%
% 4.  It computes the overall and per-segment IRI values based on this
%     simulated velocity, following the standard IRI definition.
%
% 5.  The calculated IRI values are classified and saved to new .mat and
%     .csv files in this script's own output folder.
%
% INPUTS:
% - '04_Spatial Mapping/recon_spatial.mat': The final reconstructed profile.
%
% OUTPUTS:
% - '11_IRI Reconstructed Profile/iri_recon_profile.mat': Contains the
%   calculated IRI values for the reconstructed profile.
% - '11_IRI Reconstructed Profile/iri_recon_segments.csv': A table of
%   per-segment IRI values.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ============================== HYPERPARAMETERS ==============================
% (1) Input discovery (cfg): prefer new location, fallback to legacy
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};
sp_candidates  = { ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'SpatialMapping', 'recon_spatial.mat'), ...
};

% (2) This script’s output folder
out_root   = fullfile('00_Outputs', '05_IRI_Analysis', 'ReconstructedProfile');
fig_dir    = fullfile(out_root,'figs');

% (3) Choose reconstructed series
%     'filtered' -> sp.zr_x_filt (recommended)
%     'raw'      -> sp.zr_x_raw
recon_series = 'filtered';     % 'filtered' | 'raw'

% (4) Golden Car parameters (canonical)
GC.ms_corner = 250;  % kg
GC.mu = 35;          % kg
GC.ks = 63.3e3;      % N/m
GC.cs = 6.53e3;      % Ns/m
GC.kt = 653e3;       % N/m
GC.ct = 0;           % Ns/m  (IRI standard ~0)

% (5) IRI simulation speed & numerics
iri_speed_kmh     = 80;        % km/h (IRI is defined at 80 km/h)
target_fs_time    = [];        % Hz; [] -> use native mapping t=(x-x0)/Viri
warmup_length_m   = 10;         % ignore first meters (damp IC transients)

% (6) Segmentation for reporting
segment_length_m  = 10;        % meters per IRI segment (typical 10/20/50/100 m)

% (7) Units to report and plot
%   'm_per_km' -> report in m/km only
%   'mm_per_m' -> report in mm/m only
%   'both'     -> save both in MAT/CSV; plots use m/km
iri_units_report  = 'm_per_km';     % 'm_per_km' | 'mm_per_m' | 'both'

% (8) Plots & export
make_plots  = true;
export_png  = true;

%% ------------------------------- Load inputs --------------------------------
% Pick first available cfg
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file')==2, cfg_path = cfg_candidates{k}; break; end
end
assert(~isempty(cfg_path), 'Missing recon_cfg.mat in expected locations.');

% Pick first available sp
sp_path = '';
for k = 1:numel(sp_candidates)
    if exist(sp_candidates{k},'file')==2, sp_path = sp_candidates{k}; break; end
end
assert(~isempty(sp_path), 'Missing recon_spatial.mat in expected locations.');

C  = load(cfg_path); cfg = C.cfg;
SP = load(sp_path);  sp  = SP.sp;

if ~exist(out_root,'dir'), mkdir(out_root); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Pick reconstruction series
x_rec = sp.x(:);   % distance [m]
switch lower(recon_series)
    case 'filtered', zr_rec = sp.zr_x_filt(:);
    case 'raw',      zr_rec = sp.zr_x_raw(:);
    otherwise, error('Unknown recon_series="%s"', recon_series);
end

% Ensure strictly increasing x
[x_rec, uniq_idx] = unique(x_rec, 'stable');
zr_rec = zr_rec(uniq_idx);

% IRI timing at constant speed
Viri   = iri_speed_kmh / 3.6;               % m/s
t_rec  = (x_rec - x_rec(1)) / Viri;         % s

% Optional resample in time
if ~isempty(target_fs_time)
    t_end = t_rec(end);
    t     = (0:1/target_fs_time:t_end).';
    zr    = interp1(t_rec, zr_rec, t, 'linear', 'extrap');
else
    t  = t_rec;
    zr = zr_rec;
end

% Road time derivative
% --- ADDED: 250mm Moving Average Filter (Simulate Tire Enveloping) ---
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
sys   = ss(A,B,Cout,Dout);                % continuous LTI
U     = [zr, dzr_dt];                     % inputs [zr; dzr/dt]
X     = lsim(sys, U, t, zeros(4,1));      % states: [ys, dys, yus, dyus]

ys   = X(:,1);
dys  = X(:,2);
yus  = X(:,3);
dyus = X(:,4);

% Relative suspension velocity
v_rel = dys - dyus;

% Distance aligned with t (IRI speed)
x = t * Viri + x_rec(1);

% Warmup mask and total length for normalization
m_warm = x >= (x(1) + warmup_length_m);
L_total = max(x(m_warm)) - min(x(m_warm));
L_total = max(L_total, eps);

% Overall IRI [m/km] (base unit for classification)
IRI_total_m_per_km = (1000 / L_total) * trapz(t(m_warm), abs(v_rel(m_warm)));

%% -------------------- Segment-wise IRI on fixed-length bins -----------------
segs = make_segments(x, segment_length_m);
nSeg = size(segs,1);
IRI_seg_m_per_km = nan(nSeg,1);
for k = 1:nSeg
    xs = segs(k,1); xe = segs(k,2);
    m  = (x >= xs) & (x < xe);
    if any(m)
        Lk = xe - xs; Lk = max(Lk, eps);
        IRI_seg_m_per_km(k) = (1000 / Lk) * trapz(t(m), abs(v_rel(m)));
    end
end
seg_centers = mean(segs,2);

%% ------------------------- Units & Qualification ----------------------------
% Convert to requested display units (note: 1 m/km == 1 mm/m numerically)
[overall_val, overall_units_lbl] = iri_convert_units(IRI_total_m_per_km, iri_units_report);
seg_vals = iri_convert_units(IRI_seg_m_per_km, iri_units_report);

% Classify overall and segments using base m/km thresholds
[overall_class, overall_color] = iri_classify_single(IRI_total_m_per_km);
[seg_class, seg_color, seg_class_idx] = iri_classify(IRI_seg_m_per_km);

% Also keep both-unit versions for convenience
overall_mm_per_m = IRI_total_m_per_km;   % numeric same
seg_mm_per_m     = IRI_seg_m_per_km;

%% ------------------------------- Save outputs -------------------------------
iri_rec = struct();
iri_rec.info.model        = 'Golden Car (quarter-car)';
iri_rec.info.params       = GC;
iri_rec.info.speed_kmh    = iri_speed_kmh;
iri_rec.info.warmup_m     = warmup_length_m;
iri_rec.info.seg_len_m    = segment_length_m;
iri_rec.info.recon_series = recon_series;
iri_rec.info.units_mode   = iri_units_report;

iri_rec.series.x_m        = x;
iri_rec.series.zr_m       = zr;
iri_rec.series.v_rel_mps  = v_rel;
iri_rec.series.t_s        = t;

iri_rec.overall.IRI_m_per_km  = IRI_total_m_per_km;
iri_rec.overall.IRI_mm_per_m  = overall_mm_per_m;
iri_rec.overall.value_display = overall_val;
iri_rec.overall.units_display = overall_units_lbl;
iri_rec.overall.class_name    = overall_class;
iri_rec.overall.color_rgb     = overall_color;

iri_rec.segments.edges_m        = segs;
iri_rec.segments.centers_m      = seg_centers;
iri_rec.segments.IRI_m_per_km   = IRI_seg_m_per_km;
iri_rec.segments.IRI_mm_per_m   = seg_mm_per_m;
iri_rec.segments.value_display  = seg_vals;
iri_rec.segments.class_name     = seg_class;
iri_rec.segments.class_index    = seg_class_idx;
iri_rec.segments.color_rgb      = seg_color;

save(fullfile(out_root,'iri_recon_profile.mat'), 'iri_rec','-v7.3');
fprintf('Saved: %s  (overall IRI = %.3f m/km, class = %s, series=%s)\n', ...
    fullfile(out_root,'iri_recon_profile.mat'), IRI_total_m_per_km, overall_class, recon_series);

% CSV for segments
Tseg = table((1:nSeg).', segs(:,1), segs(:,2), (segs(:,2)-segs(:,1)), ...
             IRI_seg_m_per_km, seg_mm_per_m, seg_vals(:), seg_class(:), seg_class_idx(:), ...
    'VariableNames', {'segment_idx','x_start_m','x_end_m','length_m', ...
                      'iri_m_per_km','iri_mm_per_m','iri_value_display','iri_class','class_idx'});
csv_path = fullfile(out_root,'iri_recon_segments.csv');
writetable(Tseg, csv_path);
fprintf('Saved: %s\n', csv_path);

%% --------------------------------- Plots -----------------------------------
if make_plots
    % Fig 1 — Reconstructed profile with segment grids
    f1 = mkfig(sprintf('IRI (recon:%s) — Profile & segments', recon_series), [100 70 1100 540]);
    plot(x, zr, 'b', 'LineWidth', 1.2); grid on; box on; formatAxes();
    title(sprintf('Reconstructed road profile  z_r(x)  [%s]', recon_series), 'Interpreter','none');
    xlabel('Distance (m)'); ylabel('Elevation (m)');
    yl = ylim; hold on;
    for k=1:nSeg
        plot([segs(k,1) segs(k,1)], yl, ':', 'Color',[0.7 0.7 0.7]);
    end
    if export_png, exportgraphics(f1, fullfile(fig_dir,'27_iri_recon_profile.png'), 'Resolution', 180); end

    % Fig 2 — Segment IRI bars with background bands (qualified)
    f2 = mkfig('IRI (recon) — Segment bars (qualified)', [130 90 1100 560]);
    ax = gca; hold(ax,'on'); box(ax,'on'); formatAxes(ax);

    b = bar(ax, seg_centers, seg_vals, 1.0, 'EdgeColor','none');
    b.FaceColor = 'flat';
    b.CData     = seg_color;

    yline(ax, overall_val, 'k--', ...
        sprintf('overall = %.3f %s (%s)', overall_val, overall_units_lbl, overall_class), ...
        'LabelHorizontalAlignment','left');

    % Draw bands after plotting; uses current ylim, does NOT change y-axis.
    draw_iri_bands_auto(ax, overall_units_lbl, 0.12);

    title(ax, sprintf('Segment IRI (%d m) — qualified by background bands', segment_length_m));
    xlabel(ax, 'Segment center (m)'); ylabel(ax, sprintf('IRI (%s)', overall_units_lbl));
    if export_png, exportgraphics(f2, fullfile(fig_dir,'28_iri_recon_segments_banded.png'), 'Resolution', 180); end

    % Fig 3 — Running IRI vs distance (display units)
    f3 = mkfig('IRI (recon) — Running IRI', [160 110 980 520]);
    running_mpkm = running_iri_per_km(t, x, v_rel, warmup_length_m);
    running_disp = iri_convert_units(running_mpkm, iri_units_report);
    plot(x, running_disp, 'm', 'LineWidth', 1.2); grid on; box on; formatAxes();
    title('Running IRI (from start)'); xlabel('Distance (m)'); ylabel(sprintf('IRI (%s)', overall_units_lbl));
    yline( overall_val, 'k--', 'final overall');
    if export_png, exportgraphics(f3, fullfile(fig_dir,'29_iri_recon_running.png'), 'Resolution', 180); end

    % Fig 4 — Relative suspension velocity trace (warmup highlighted)
    f4 = mkfig('IRI (recon) — Relative suspension velocity', [190 130 1100 520]);
    plot(x, v_rel, 'Color',[0.25 0.45 0.85], 'LineWidth', 1.0); grid on; box on; formatAxes(); hold on;
    title('Relative suspension velocity  v_{rel}(x) = z_ṡ - z_u̇');
    xlabel('Distance (m)'); ylabel('Velocity (m/s)');
    xline(x(1)+warmup_length_m, 'r--', sprintf('warmup = %dm', warmup_length_m));
    if export_png, exportgraphics(f4, fullfile(fig_dir,'30_iri_recon_vrel.png'), 'Resolution', 180); end
end


%% =============================== FUNCTIONS ===============================
function [A,B] = quarter_AB(p)
% Quarter-car (Golden Car style) continuous-time state-space:
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
    if n<2, return; end
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
    if edges(end) < xmax, edges = [edges xmax]; end
    segs = [edges(1:end-1).' edges(2:end).'];
end

function iri_running = running_iri_per_km(t, x, v_rel, warmup_m)
% Compute cumulative IRI along the path; ignore first warmup_m meters.
    iri_running = nan(size(x));
    x0 = x(1) + warmup_m;     % use function arg, not outer var
    m0 = x >= x0;
    if ~any(m0)
        iri_running(:) = 0; return;
    end
    idx0 = find(m0,1,'first');
    cumI = cumtrapz(t(idx0:end), abs(v_rel(idx0:end)));
    Ls   = x(idx0:end) - x(idx0);
    Ls(Ls<=0) = eps;
    iri_running(idx0:end) = 1000 * (cumI ./ Ls);  % m/km
    iri_running(1:idx0-1) = NaN;
end

% ---------- Qualification helpers (same thresholds for m/km and mm/m) ----------
function [val_disp, units_lbl] = iri_convert_units(val_m_per_km, mode)
    switch lower(mode)
        case 'm_per_km'
            val_disp = val_m_per_km; units_lbl = 'm/km';
        case 'mm_per_m'
            val_disp = val_m_per_km; units_lbl = 'mm/m';
        case 'both'
            val_disp = val_m_per_km; units_lbl = 'm/km';
        otherwise
            error('Unknown iri_units_report: %s', mode);
    end
end

function [name, rgb] = iri_classify_single(iri_mpkm)
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
    name = names{idx}; rgb = colors(idx,:);
end

function [names, colors, idx] = iri_classify(iri_mpkm_vec)
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
    colors = colorsLUT(idx,:);
end

function [names, colors, edges] = iri_classes_def()
    names  = {'Very Good','Good','Fair','Poor','Very Poor'};
    colors = [  0.00 0.60 0.05;   % green
                0.70 0.82 0.20;   % yellow-green
                0.95 0.70 0.20;   % orange
                0.90 0.00 0.00;   % red
                0.60 0.10 0.10];  % dark red
    edges  = [0 2 4 6 10];        % m/km thresholds
end

function hBands = draw_iri_bands_auto(ax, units_lbl, alpha_bg)
% Draw semi-transparent IRI class bands WITHOUT changing y-limits.
% Place this AFTER you plot bars/lines so the axis has auto-scaled.
% Returns array of Patch handles (already pushed behind other content).

    if nargin < 3 || isempty(alpha_bg), alpha_bg = 0.12; end

    [~, baseColors, edges] = iri_classes_def();  % edges = [0 2 4 6 10] (m/km)
    colors = 1 - (1 - baseColors) * 0.35;        % lighten for backgrounds

    % Use current axis limits (do NOT modify them)
    xl = xlim(ax);
    yl = ylim(ax);

    % Upper edge of each band: [2 4 6 10  (current top)]
    uppers = [edges(2:end) yl(2)];
    lowers = edges;

    % Draw only the part of each band that actually lies inside current ylim
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

    % Threshold lines (only if they fall within current ylim)
    for e = edges(2:end)
        if e > yl(1) && e < yl(2)
            yline(ax, e, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75, 'HandleVisibility','off');
        end
    end

    % Compact labels at the right edge (inside axes)
    DX = diff(xl);
    x_text = xl(2) - 0.01*DX;
    labelTxt = {sprintf('< 2.0 %s', units_lbl), sprintf('2.0–4.0 %s', units_lbl), ...
                sprintf('4.0–6.0 %s', units_lbl), sprintf('6.0–10.0 %s', units_lbl), ...
                sprintf('> 10.0 %s', units_lbl)};
    for i = 1:5
        y0 = max(yl(1), lowers(i));
        y1 = min(yl(2), uppers(i));
        if y1 > y0
            y_mid = (y0 + y1)/2;
            text(ax, x_text, y_mid, labelTxt{i}, ...
                'HorizontalAlignment','right','VerticalAlignment','middle', ...
                'FontSize',8,'Color',[0.25 0.25 0.25], 'HitTest','off');
        end
    end

    % Keep valid bands behind everything else
    hb = hBands(isgraphics(hBands));
    if ~isempty(hb), uistack(hb, 'bottom'); end
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