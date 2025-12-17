
%% 12) Compare IRI of True vs. Reconstructed Profiles
% -------------------------------------------------------------------------
% PURPOSE:
% This script directly compares the International Roughness Index (IRI)
% calculated from the ground truth profile against the IRI calculated from
% the reconstructed profile. It provides the ultimate validation of the
% reconstruction's performance in terms of this key industry metric.
%
% CORE LOGIC (How it works):
% 1.  Loads the two main IRI result files: 'iri_original_profile.mat' (from
%     the true profile) and 'iri_recon_profile.mat' (from the reconstructed
%     profile).
%
% 2.  To ensure a fair comparison, it identifies the common overlapping
%     distance range between the two profiles and establishes a common grid
%     of road segments.
%
% 3.  It re-calculates the per-segment IRI values for BOTH the true and
%     reconstructed profiles on this common grid, ensuring an apples-to-
%     apples comparison for each segment.
%
% 4.  It computes the error (absolute and relative) between the true and
%     reconstructed IRI, both for the overall value and for each segment.
%
% 5.  Saves all comparison metrics to .mat and .csv files.
%
% 6.  Generates several comparison plots, such as a grouped bar chart
%     showing the true vs. reconstructed IRI for each segment, and a
%     scatter plot to show their correlation.
%
% INPUTS:
% - '10_IRI True Profile/iri_original_profile.mat'
% - '11_IRI Reconstructed Profile/iri_recon_profile.mat'
%
% OUTPUTS:
% - '12_IRI Compare/iri_compare.mat': Contains all comparison metrics.
% - '12_IRI Compare/iri_compare_segments.csv': A detailed per-segment table
%   of true IRI, reconstructed IRI, and the error.
% - Figures visualizing the comparison.
% -------------------------------------------------------------------------

close all; clc;

%% ============================== HYPERPARAMETERS ==============================
% (1) Input discovery (prefer new folders, fallback to legacy)
if SIM_SOURCE == "LTPP"
    orig_candidates  = {fullfile('00_Outputs', '05_IRI_Analysis', 'TrueProfile', 'iri_true_profile.mat')};
else 
    orig_candidates  = {fullfile('00_Outputs', '05_IRI_Analysis', 'TrueProfile', 'iri_original_profile.mat')};
end

recon_candidates = { ...
    fullfile('00_Outputs', '05_IRI_Analysis', 'ReconstructedProfile', 'iri_recon_profile.mat'), ...
};

% (2) Outputs
out_root   = fullfile('00_Outputs', '05_IRI_Analysis', 'Comparison');
fig_dir    = fullfile(out_root,'figs');

% (3) Segmenting & warmup
segment_length_m    = 10;   % segment size for comparison
common_warmup_m     = 0;   % ignore first X meters from each series' start

% (4) Units for reporting/plots
%   'm_per_km' -> only m/km
%   'mm_per_m' -> only mm/m
%   'both'     -> save both; plots use m/km label
iri_units_report    = 'm_per_km';  % 'm_per_km' | 'mm_per_m' | 'both'

% (5) Plots & export
make_plots = true;
export_png = true;

%% ------------------------------ Load inputs ---------------------------------
orig_path  = first_existing(orig_candidates);
recon_path = first_existing(recon_candidates);
assert(~isempty(orig_path),  'iri_original_profile.mat not found in expected locations.');
assert(~isempty(recon_path), 'iri_recon_profile.mat not found in expected locations.');

Strue = load(orig_path);
Sreco = load(recon_path);

% Accept either computed 'iri' (series) or imported 'IRI' (segments from ProVAL)
if isfield(Strue,'iri') && isfield(Strue.iri,'series')
    iri = Strue.iri;
    true_mode = "series";
elseif isfield(Strue,'IRI') && isfield(Strue.IRI,'segments')
    iri_imp = Strue.IRI;
    true_mode = "segments";
else
    error('Could not find a valid true-IRI struct in %s (need ''iri.series'' or ''IRI.segments'').', orig_path);
end

iri_rec = Sreco.iri_rec;

if ~exist(out_root,'dir')
    mkdir(out_root);
end
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Sanity
assert(isfield(iri_rec,'series') && isfield(iri_rec.series,'x_m') && isfield(iri_rec.series,'v_rel_mps'), 'Bad iri_rec struct.');

if strcmp(true_mode,"series")
    assert(isfield(iri,'series') && isfield(iri.series,'x_m') && isfield(iri.series,'v_rel_mps'), 'Bad iri struct (series).');
else
    Ttrue = iri_imp.segments;
    req = all(ismember({'Start_m','End_m','Length_m','IRI_m_per_km'}, Ttrue.Properties.VariableNames));
    assert(req, 'Bad imported IRI segments table.');
end

% Pull TRUE data
if strcmp(true_mode,"series")
    xT = iri.series.x_m(:);
    tT = iri.series.t_s(:);
    vT = iri.series.v_rel_mps(:);
else
    xT = []; tT = []; vT = [];
    % Imported segment edges and values
    Ttrue  = iri_imp.segments;
    sEdges = [Ttrue.Start_m, Ttrue.End_m];
    sIRI   = Ttrue.IRI_m_per_km;
end

% Pull series (RECON)
xR = iri_rec.series.x_m(:);
tR = iri_rec.series.t_s(:);
vR = iri_rec.series.v_rel_mps(:);

% Overlap region
if strcmp(true_mode,"series")
    [xmin, xmax] = safe_overlap_range([xT(1) xT(end)], [xR(1) xR(end)]);
else
    rngTrue = [min(sEdges(:,1)) max(sEdges(:,2))];
    [xmin, xmax] = safe_overlap_range(rngTrue, [xR(1) xR(end)]);
end
assert(xmax > xmin, 'No spatial overlap between original and reconstructed paths.');

% Warmup
if strcmp(true_mode,"series")
    x0T = xT(1) + common_warmup_m;
else
    x0T = xmin + common_warmup_m;
end
x0R = xR(1) + common_warmup_m;

% Masks for recon (and true if series)
mR = (xR >= max(xmin, x0R)) & (xR <= xmax);
assert(any(mR), 'Warmup + overlap left no recon data to compare.');
if strcmp(true_mode,"series")
    mT = (xT >= max(xmin, x0T)) & (xT <= xmax);
    assert(any(mT), 'Warmup + overlap left no true data to compare.');
end

% Build segment edges on the qualified window
keep = [];
if strcmp(true_mode,"series")
    xStart = max(xmin, x0T);
    segs = make_segments_common(xStart, xmax, segment_length_m);
else
    xStart = max(xmin, x0T);
    E = sEdges;
    keep = (E(:,2) > xStart) & (E(:,1) < xmax);
    E = E(keep,:);
    E(:,1) = max(E(:,1), xStart);
    E(:,2) = min(E(:,2), xmax);
    segs = E;
end
seg_centers = mean(segs,2);

%% -------------------- Overall IRI (on the common overlap) -------------------
if strcmp(true_mode,"series")
    L_true = max(xT(mT)) - min(xT(mT));  L_true = max(L_true, eps);
    L_reco = max(xR(mR)) - min(xR(mR));  L_reco = max(L_reco, eps);

    IRI_true_overall = 1000 * trapz(tT(mT), abs(vT(mT))) / L_true;   % m/km
    IRI_reco_overall = 1000 * trapz(tR(mR), abs(vR(mR))) / L_reco;   % m/km

    overall_abs_err = IRI_reco_overall - IRI_true_overall;
    overall_rel_err = overall_abs_err / max(IRI_true_overall, eps); % signed fraction

    [true_overall_disp, units_lbl] = iri_convert_units(IRI_true_overall, iri_units_report);
    [reco_overall_disp, ~] = iri_convert_units(IRI_reco_overall, iri_units_report);

    [true_overall_class,  true_overall_color] = iri_classify_single(IRI_true_overall);
    [reco_overall_class,  reco_overall_color] = iri_classify_single(IRI_reco_overall);
else
    % Defer overall computation to the segment stage (uses imported segments)
    IRI_true_overall = NaN; IRI_reco_overall = NaN;
    overall_abs_err = NaN;  overall_rel_err = NaN;
    units_lbl = 'm/km';
    true_overall_disp = NaN; reco_overall_disp = NaN;
    true_overall_class = ''; true_overall_color = [0 0 0];
    reco_overall_class = ''; reco_overall_color = [0 0 0];
end

%% ---------------- Segment-wise IRI (rebuilt on common segments) -------------
if strcmp(true_mode,"series")
    IRI_true_seg = segment_iri_from_series(segs, xT, tT, vT);   % m/km
else
    % Use imported per-segment IRIs aligned to kept edges
    IRI_true_seg = sIRI(keep);
end
IRI_reco_seg = segment_iri_from_series(segs, xR, tR, vR);   % m/km

% If overall was deferred (segments mode), compute it now from segs
if isnan(IRI_true_overall)
    len_seg = segs(:,2) - segs(:,1);
    IRI_true_overall = sum(IRI_true_seg .* len_seg) / sum(len_seg);
    IRI_reco_overall = sum(IRI_reco_seg .* len_seg) / sum(len_seg);

    overall_abs_err = IRI_reco_overall - IRI_true_overall;
    overall_rel_err = overall_abs_err / max(IRI_true_overall, eps);

    [true_overall_disp, units_lbl] = iri_convert_units(IRI_true_overall, iri_units_report);
    [reco_overall_disp, ~] = iri_convert_units(IRI_reco_overall, iri_units_report);

    [true_overall_class,  true_overall_color] = iri_classify_single(IRI_true_overall);
    [reco_overall_class,  reco_overall_color] = iri_classify_single(IRI_reco_overall);
end

% Errors
seg_abs_err = IRI_reco_seg - IRI_true_seg;         % m/km
seg_rel_err = seg_abs_err ./ max(IRI_true_seg, eps);

% Classify segments
[seg_class_true, seg_color_true, seg_idx_true] = iri_classify(IRI_true_seg);
[seg_class_reco, seg_color_reco, seg_idx_reco] = iri_classify(IRI_reco_seg);

% Units for display
seg_true_disp = iri_convert_units(IRI_true_seg, iri_units_report);
seg_reco_disp = iri_convert_units(IRI_reco_seg, iri_units_report);

% Goodness
[r_seg, a1, a0] = simple_linreg(IRI_true_seg, IRI_reco_seg);   % slope, intercept, r

%% ------------------------------ Save outputs --------------------------------
C = struct();
C.settings.segment_length_m = segment_length_m;
C.settings.common_warmup_m = common_warmup_m;
C.settings.units_mode = iri_units_report;

C.overlap.xmin = xmin;
C.overlap.xmax = xmax;

C.overall.true_m_per_km = IRI_true_overall;
C.overall.reco_m_per_km = IRI_reco_overall;
C.overall.abs_err_m_per_km = overall_abs_err;
C.overall.rel_err          = overall_rel_err;

C.overall.true_display_value = true_overall_disp;
C.overall.reco_display_value = reco_overall_disp;
C.overall.units_display      = units_lbl;
C.overall.true_class_name    = true_overall_class;
C.overall.true_color_rgb     = true_overall_color;
C.overall.reco_class_name    = reco_overall_class;
C.overall.reco_color_rgb     = reco_overall_color;

C.segments.edges_m          = segs;
C.segments.centers_m        = seg_centers;
C.segments.true_m_per_km    = IRI_true_seg;
C.segments.reco_m_per_km    = IRI_reco_seg;
C.segments.abs_err_m_per_km = seg_abs_err;
C.segments.rel_err          = seg_rel_err;

C.segments.true_display     = seg_true_disp;
C.segments.reco_display     = seg_reco_disp;
C.segments.units_display    = units_lbl;
C.segments.true_class       = seg_class_true;
C.segments.true_class_idx   = seg_idx_true;
C.segments.true_color_rgb   = seg_color_true;
C.segments.reco_class       = seg_class_reco;
C.segments.reco_class_idx   = seg_idx_reco;
C.segments.reco_color_rgb   = seg_color_reco;

save(fullfile(out_root,'iri_compare.mat'), 'C','-v7.3');
fprintf('Saved: %s\n', fullfile(out_root,'iri_compare.mat'));

% CSV (segment-wise)
T = table((1:size(segs,1)).', segs(:,1), segs(:,2), seg_centers, ...
          IRI_true_seg, IRI_reco_seg, seg_abs_err, seg_rel_err, ...
          seg_true_disp, seg_reco_disp, ...
          seg_class_true(:), seg_idx_true(:), seg_class_reco(:), seg_idx_reco(:), ...
    'VariableNames', {'segment_idx','x_start_m','x_end_m','x_center_m', ...
                      'iri_true_m_per_km','iri_reco_m_per_km', ...
                      'abs_err_m_per_km','rel_err', ...
                      'iri_true_display','iri_reco_display', ...
                      'class_true','class_idx_true','class_reco','class_idx_reco'});
csv_path = fullfile(out_root,'iri_compare_segments.csv');
writetable(T, csv_path);
fprintf('Saved: %s\n', csv_path);

%% --------------------------------- Plots ------------------------------------
if make_plots
    % Fig 1 — Segment bars (grouped) with background class bands
    f1 = mkfig('IRI compare — segments (qualified)', [100 70 1200 720]);
    ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); formatAxes(ax);

    B = bar(ax, seg_centers, [seg_true_disp(:) seg_reco_disp(:)], 'grouped');
    pastel_true = 1 - (1 - seg_color_true) * 0.55;
    B(1).FaceColor='flat'; B(1).CData=pastel_true;   B(1).EdgeColor='none'; B(1).FaceAlpha=0.9;
    B(2).FaceColor='flat'; B(2).CData=seg_color_reco; B(2).EdgeColor='none'; B(2).FaceAlpha=1.0;

    % Overall lines (already in display units)
    yline(ax, true_overall_disp, '--', sprintf('truth = %.3f %s (%s)', ...
        true_overall_disp, units_lbl, true_overall_class), ...
        'Color',[0.3 0.3 0.3], 'LabelHorizontalAlignment','left');
    yline(ax, reco_overall_disp, '--', sprintf('recon = %.3f %s (%s)', ...
        reco_overall_disp, units_lbl, reco_overall_class), ...
        'Color',[0.1 0.35 0.9], 'LabelHorizontalAlignment','left');

    % Draw bands AFTER bars so y-limits remain fully auto
    draw_iri_bands_auto(ax, units_lbl, 0.12);

    title(ax, sprintf('Segment IRI (L_{seg} = %d m) — qualified background', segment_length_m));
    xlabel(ax, 'Segment center (m)'); ylabel(ax, sprintf('IRI (%s)', units_lbl));
    legend(ax, {'Truth','Recon'}, 'Location','northwest');

    if export_png
        exportgraphics(f1, fullfile(fig_dir,'31_iri_segment_bars.png'), 'Resolution', 180);
    end

    % Fig 2 — Scatter (Recon vs Truth) colored by TRUTH class (always m/km)
    f2 = mkfig('IRI compare — scatter (qualified)', [140 100 900 700]);
    hold on; grid on; box on; formatAxes();
    [~, colors] = iri_classes_def();
    for k = 1:5
        m = seg_idx_true == k;
        scatter(IRI_true_seg(m), IRI_reco_seg(m), 28, colors(k,:), ...
            'filled', 'MarkerEdgeColor','k','MarkerFaceAlpha',0.9);
    end
    xr = [min(IRI_true_seg), max(IRI_true_seg)];
    xr(1) = xr(1) - 0.02*range(xr); xr(2) = xr(2) + 0.02*range(xr);
    plot(xr, xr, 'k--', 'LineWidth', 1.0);
    plot(xr, a1*xr + a0, 'r-', 'LineWidth', 1.2);
    title(sprintf('Segment IRI: recon vs truth   (r=%.3f, y=%.3fx%+ .3f)', r_seg, a1, a0));
    xlabel('Truth IRI (m/km)'); ylabel('Recon IRI (m/km)');
    if export_png
        exportgraphics(f2, fullfile(fig_dir,'32_iri_scatter_segments.png'), 'Resolution', 180);
    end

    % Fig 3 — Running IRI (only when TRUE has time-series)
    if strcmp(true_mode,"series")
        f3 = mkfig('IRI compare — running IRI', [170 120 980 540]);
        runT = running_iri_on_range(xT, tT, vT, max(xmin, x0T));
        runR = running_iri_on_range(xR, tR, vR, max(xmin, x0R));
        plot(runT.x, iri_convert_units(runT.iri_m_per_km, iri_units_report), 'Color',[0.2 0.2 0.2], 'LineWidth', 1.3); hold on;
        plot(runR.x, iri_convert_units(runR.iri_m_per_km, iri_units_report), 'Color',[0.1 0.35 0.9], 'LineWidth', 1.3);
        grid on; box on; formatAxes();
        title('Running IRI from start of overlap (post-warmup)');
        xlabel('Distance (m)'); ylabel(sprintf('IRI (%s)', units_lbl));
        legend('Truth','Recon','Location','best');
        if export_png
            exportgraphics(f3, fullfile(fig_dir,'33_iri_running_compare.png'), 'Resolution', 180);
        end
    end

    % Fig 4 — Overall summary table
    f4 = mkfig('IRI compare — overall summary', [200 160 760 380]);
    axis off;
    lines = {
        sprintf('Overlap: [%.1f, %.1f] m', xmin, xmax)
        sprintf('Segment length: %d m    |  Warmup: %d m', segment_length_m, common_warmup_m)
        ' '
        sprintf('Overall IRI (truth):  %.3f m/km  (%s)', IRI_true_overall, true_overall_class)
        sprintf('Overall IRI (recon):  %.3f m/km  (%s)', IRI_reco_overall, reco_overall_class)
        sprintf('Abs error:            %+ .3f m/km', overall_abs_err)
        sprintf('Rel error:            %+ .2f %%', 100*overall_rel_err)
        ' '
        sprintf('Segment correlation r: %.3f', r_seg)
        sprintf('Segment fit:           y = %.3f x %+ .3f', a1, a0)
    };
    text(0.07, 0.92, lines, 'Units','normalized', 'FontName','Consolas','FontSize',11);
    if export_png
        exportgraphics(f4, fullfile(fig_dir,'34_iri_overall_table.png'), 'Resolution', 180);
    end
end

%% =============================== FUNCTIONS ===============================
function p = first_existing(candidates)
    p = '';
    for i = 1:numel(candidates)
        if exist(candidates{i},'file')==2, p = candidates{i}; return; end
    end
end

function [xmin, xmax] = safe_overlap_range(r1, r2)
    xmin = max(min(r1), min(r2));
    xmax = min(max(r1), max(r2));
end

function segs = make_segments_common(xmin, xmax, Lseg)
% Build [start,end] edges over [xmin,xmax] with fixed length Lseg.
    x0 = xmin;
    edges = x0 : Lseg : (xmax + 1e-9);
    if edges(end) < xmax, edges = [edges xmax]; end
    segs = [edges(1:end-1).' edges(2:end).'];
end

function I = segment_iri_from_series(segs, x, t, vrel)
% Per-segment IRI = (1000/Lseg) * ∫ |v_rel| dt, with samples inside segment.
    nSeg = size(segs,1);
    I = nan(nSeg,1);
    for k = 1:nSeg
        xs = segs(k,1); xe = segs(k,2);
        m  = (x >= xs) & (x < xe);
        if any(m)
            Lk = xe - xs; Lk = max(Lk, eps);
            I(k) = (1000 / Lk) * trapz(t(m), abs(vrel(m)));
        end
    end
end

function R = running_iri_on_range(x, t, vrel, xstart)
% Running IRI from xstart to end of x; returns struct with x & IRI (m/km).
    m0 = x >= xstart;
    if ~any(m0)
        R.x = x; R.iri_m_per_km = nan(size(x)); return;
    end
    x1 = x(find(m0,1,'first'));
    idx0 = find(m0,1,'first');
    cumI = cumtrapz(t(idx0:end), abs(vrel(idx0:end)));
    Ls   = x(idx0:end) - x1;
    Ls(Ls <= 0) = eps;
    R.x = x(idx0:end);
    R.iri_m_per_km = 1000 * (cumI ./ Ls);
end

function [r, a1, a0] = simple_linreg(x, y)
% Pearson r and least-squares fit y = a1*x + a0.
    x = x(:); y = y(:);
    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    if numel(x) < 3
        r = NaN; a1 = NaN; a0 = NaN; return;
    end
    X = [x ones(size(x))];
    beta = X \ y;
    a1 = beta(1); a0 = beta(2);
    C = corrcoef(x, y);
    r = C(1,2);
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

    % Optional compact labels at the right edge (inside axes)
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

    % Keep bands behind everything else (only valid handles)
    hb = hBands(isgraphics(hBands));
    if ~isempty(hb), uistack(hb, 'bottom'); end
    ax.Layer = 'top';
end