
%% 8) Generate Visual Comparison Plots vs. Ground Truth
% -------------------------------------------------------------------------
% PURPOSE:
% This script provides a visual comparison between the final reconstructed
% road profile and the original ground truth. It generates a series of
% plots to intuitively illustrate the accuracy and error characteristics of
% the reconstruction.
%
% CORE LOGIC (How it works):
% 1.  Loads the reconstructed spatial profile ('recon_spatial.mat') and the
%     configuration file ('recon_cfg.mat'), which contains the ground truth.
%
% 2.  Generates an overlay plot showing the reconstructed profile drawn
%     directly on top of the ground truth profile, along with a separate
%     plot of the error (the difference between the two) at each point.
%
% 3.  Creates a histogram and an empirical Cumulative Distribution Function
%     (CDF) of the error signal to visualize the distribution of errors.
%
% 4.  Generates a comparison of the Power Spectral Density (PSD) of both
%     profiles, showing how well the frequency content of the road was
%     reconstructed.
%
% 5.  Automatically identifies areas with significant road features (like
%     bumps or potholes) and creates zoomed-in overlay plots for a detailed
%     look at performance in these critical areas.
%
% 6.  Creates a scatter plot of reconstructed elevation vs. true elevation
%     to visualize correlation and any systematic bias.
%
% 7.  Saves all generated figures as image files.
%
% INPUTS:
% - '02_Reconstruction Inputs/recon_cfg.mat': For the ground truth profile.
% - '04_Spatial Mapping/recon_spatial.mat': For the reconstructed profile.
%
% OUTPUTS:
% - Image files (e.g., .png) saved to the '08_Comparison Plots/figs' folder.
% -------------------------------------------------------------------------

close all; clear; clc;

%% ========================== HYPERPARAMETERS ==========================
% (1) Input discovery (new locations first, then legacy)
cfg_candidates = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};
sp_candidates  = { ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'SpatialMapping', 'recon_spatial.mat'), ...
};

% (2) This script’s output folder (per-script organization)
save_folder = fullfile('00_Outputs', '04_Validation', 'ComparisonPlots');
fig_dir     = fullfile(save_folder,'figs');

% (3) Which reconstructed series to plot?
%     'filtered' -> use zr_x_filt (recommended; paper low-pass policy)
%     'raw'      -> use zr_x_raw  (pre-filter)
recon_series = 'filtered';   % 'filtered' | 'raw'

% (4) Error zoom windows (auto-detected around strong features)
n_zoom_windows        = 3;     % how many zoom panels
zoom_window_length_m  = 15.0;  % window length in meters for each zoom
min_peak_distance_m   = 30.0;  % minimum separation between zoom centers

% (5) Spatial PSD settings
welch_nfft_points  = [];       % [] => auto; else integer
welch_window_frac  = 0.50;     % window length as fraction of N
welch_overlap_frac = 0.50;     % 50% overlap
remove_dc_for_psd  = true;     % remove mean before PSD

% (6) Plots & export
make_plots = true;
export_png = true;
% =====================================================================

%% ------------------------- Locate & load inputs -------------------------
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file')==2, cfg_path = cfg_candidates{k}; break; end
end
assert(~isempty(cfg_path), 'Missing recon_cfg.mat in expected locations.');

sp_path = '';
for k = 1:numel(sp_candidates)
    if exist(sp_candidates{k},'file')==2, sp_path = sp_candidates{k}; break; end
end
assert(~isempty(sp_path), 'Missing recon_spatial.mat in expected locations.');

% Make output folders
if ~exist(save_folder,'dir'), mkdir(save_folder); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

C  = load(cfg_path); cfg = C.cfg;
SP = load(sp_path);  sp  = SP.sp;

% Ground truth (required)
assert(isfield(cfg,'ground') && isfield(cfg.ground,'x') && ~isempty(cfg.ground.x), ...
    'Ground-truth profile not found in cfg.ground.x. Provide ground truth to plot.');
x_true  = cfg.ground.x(:);
zr_true = cfg.ground.zr(:);

% Reconstructed series
x_rec   = sp.x(:);
switch lower(recon_series)
    case 'filtered', zr_rec = sp.zr_x_filt(:);
    case 'raw',      zr_rec = sp.zr_x_raw(:);
    otherwise, error('Unknown recon_series="%s"', recon_series);
end

% Common support
xmin = max(min(x_true), min(x_rec));
xmax = min(max(x_true), max(x_rec));
mT = (x_true>=xmin) & (x_true<=xmax);
mR = (x_rec >=xmin) & (x_rec <=xmax);
xT = x_true(mT); yT = zr_true(mT);
xR = x_rec(mR);  yR = zr_rec(mR);

% Interpolate reconstructed onto truth grid for direct comparison
y_rec_on_true = safe_interp(xR, yR, xT);
err = y_rec_on_true - yT;

% Spatial sampling
dx         = mean(diff(xT),'omitnan');
fs_spatial = 1 / max(dx, eps);       % samples per meter
nyq_sp     = fs_spatial / 2;         % cycles/m
nc         = NaN; if isfield(sp,'nc_spatial'), nc = sp.nc_spatial; end

%% ------------------------- PLOTS -------------------------
if make_plots
    %% Fig 1 — Overlay & error (full length)
    f1 = mkfig('Truth vs Reconstruction (full) + Error', [100 70 1150 780]);
    tlo1 = tiledlayout(f1,3,1,'TileSpacing','compact','Padding','compact');

    nexttile(tlo1,1);
    plot(xT, yT, 'k', 'LineWidth', 1.2); hold on;
    plot(xT, y_rec_on_true, 'b', 'LineWidth', 1.2);
    grid on; box on; formatAxes();
    title(sprintf('Road Elevation — Truth vs Recon (%s)', recon_series));
    xlabel('Distance (m)'); ylabel('Elevation (m)');
    legend('Truth','Recon','Location','best');

    nexttile(tlo1,2);
    plot(xT, err, 'Color',[0.25 0.45 0.85], 'LineWidth', 1.0);
    grid on; box on; formatAxes();
    title('Error e(x) = z_r^{rec}(x) - z_r^{true}(x)');
    xlabel('Distance (m)'); ylabel('Error (m)'); yline(0,'k:');

    nexttile(tlo1,3);
    % Running RMS of error (sliding window ~ 50 m or 5% length, whichever bigger)
    L = xT(end)-xT(1);
    wlen_m = max(50, 0.05*L);
    wlen_samp = max(5, round(wlen_m * fs_spatial));
    err_rms = sqrt(movmean(err.^2, wlen_samp, 'Endpoints','shrink'));
    plot(xT, err_rms, 'm', 'LineWidth', 1.1); grid on; box on; formatAxes();
    title(sprintf('Running RMS error (window ~ %.0f m)', wlen_m));
    xlabel('Distance (m)'); ylabel('RMS error (m)');

    if export_png, exportgraphics(f1, fullfile(fig_dir,'15_overlay_and_error.png'), 'Resolution', 180); end

    %% Fig 2 — Error histogram & empirical CDF
    f2 = mkfig('Error histogram & CDF', [140 90 1050 520]);

    subplot(1,2,1);
    histogram(err, 60, 'FaceColor',[0.75 0.85 1.0], 'EdgeColor',[0.2 0.2 0.6]); grid on; box on; formatAxes();
    title('Error histogram'); xlabel('Error (m)'); ylabel('Count');

    subplot(1,2,2);
    [f_cdf, x_cdf] = ecdf(err);
    plot(x_cdf, f_cdf, 'b', 'LineWidth', 1.2); grid on; box on; formatAxes();
    title('Empirical CDF of error'); xlabel('Error (m)'); ylabel('F(e)');

    if export_png, exportgraphics(f2, fullfile(fig_dir,'16_error_hist_cdf.png'), 'Resolution', 180); end

    %% Fig 3 — Spatial PSDs and PSD ratio
    f3 = mkfig('Spatial PSDs & ratio', [160 110 1100 700]);
    tlo3 = tiledlayout(f3,3,1,'TileSpacing','compact','Padding','compact');

    [Ptrue, Ntrue, welch_cfg] = welch_psd_spatial(yT, fs_spatial, ...
                        welch_nfft_points, welch_window_frac, welch_overlap_frac, remove_dc_for_psd);
    [Prec,  Nrec]             = welch_psd_spatial(y_rec_on_true, fs_spatial, ...
                        welch_nfft_points, welch_window_frac, welch_overlap_frac, remove_dc_for_psd);

    nexttile(tlo3,1);
    loglog(Ntrue, Ptrue, 'k', 'LineWidth', 1.1); hold on;
    loglog(Nrec,  Prec,  'b', 'LineWidth', 1.1); grid on; box on; formatAxes();
    title('Spatial PSD — Truth vs Recon'); xlabel('Spatial frequency (cycles/m)'); ylabel('PSD (m^2·m)');
    if isfinite(nc), add_vline(nc, '-', sprintf('LP cutoff n_c=%.3f', nc)); end
    legend('Truth','Recon','Location','southwest');

    % Common grid and PSD ratio
    n_common = linspace(0, min([Ntrue(end), Nrec(end), nyq_sp]), max(numel(Ntrue),numel(Nrec)))';
    Ptrue_c  = safe_interp(Ntrue, Ptrue, n_common);
    Prec_c   = safe_interp(Nrec,  Prec,  n_common);
    ratio_psd = Prec_c ./ max(Ptrue_c, eps);

    nexttile(tlo3,2);
    semilogx(n_common, 10*log10(Ptrue_c+eps), 'k', 'LineWidth', 1.0); hold on;
    semilogx(n_common, 10*log10(Prec_c+eps),  'b', 'LineWidth', 1.0);
    grid on; box on; formatAxes();
    title('Spatial PSD (dB)'); xlabel('Spatial frequency (cycles/m)'); ylabel('PSD (dB re m^2·m)');
    if isfinite(nc), add_vline(nc, '-', 'n_c'); end
    legend('Truth','Recon','Location','southwest');

    nexttile(tlo3,3);
    semilogx(n_common, ratio_psd, 'm', 'LineWidth', 1.1); grid on; box on; formatAxes();
    title('PSD ratio  Recon/Truth'); xlabel('Spatial frequency (cycles/m)'); ylabel('ratio');
    yline(1,'k:'); if isfinite(nc), add_vline(nc,'-','n_c'); end

    if export_png, exportgraphics(f3, fullfile(fig_dir,'17_spatial_psd_and_ratio.png'), 'Resolution', 180); end

    %% Fig 4 — Zoom windows around strong features
    f4 = mkfig('Zoomed overlays around features', [180 130 1200 860]);

    % Find prominent feature centers based on |second derivative| of truth
    y2 = second_derivative(yT, dx);
    [center_idx] = pick_zoom_centers(xT, abs(y2), n_zoom_windows, min_peak_distance_m);
    rows = n_zoom_windows; cols = 1;
    tlo4 = tiledlayout(f4, rows, cols, 'TileSpacing','compact','Padding','compact');

    for i=1:numel(center_idx)
        xi = xT(center_idx(i));
        xs = max(xT(1), xi - zoom_window_length_m/2);
        xe = min(xT(end), xi + zoom_window_length_m/2);

        m = (xT>=xs) & (xT<=xe);
        nexttile(tlo4,i);
        plot(xT(m), yT(m), 'k', 'LineWidth', 1.2); hold on;
        plot(xT(m), y_rec_on_true(m), 'b', 'LineWidth', 1.2);
        grid on; box on; formatAxes();
        title(sprintf('Zoom %d: x \\in [%.1f, %.1f] m', i, xs, xe), 'Interpreter','tex');
        xlabel('Distance (m)'); ylabel('Elevation (m)');
        legend('Truth','Recon','Location','best');
    end

    if export_png, exportgraphics(f4, fullfile(fig_dir,'18_zoom_windows.png'), 'Resolution', 180); end

    %% Fig 5 — Scatter plot with regression
    f5 = mkfig('Scatter: Truth vs Recon', [210 160 780 680]);
    scatter(yT, y_rec_on_true, 8, [0.25 0.45 0.85], 'filled'); grid on; box on; formatAxes();
    hold on;
    % Regression line
    [a1, a0, r_val] = simple_linreg(yT, y_rec_on_true);
    xr = [min(yT), max(yT)];
    yr = a1*xr + a0;
    plot(xr, yr, 'r-', 'LineWidth', 1.4);
    % y=x line
    plot(xr, xr, 'k--', 'LineWidth', 1.0);
    title('Scatter of z_r^{true} vs z_r^{rec}');
    xlabel('Truth (m)'); ylabel('Recon (m)');
    legend('samples', sprintf('fit: y=%.3fx%+ .3f (r=%.3f)', a1, a0, r_val), 'y=x', 'Location','best');

    if export_png, exportgraphics(f5, fullfile(fig_dir,'19_scatter_truth_vs_recon.png'), 'Resolution', 180); end
end

%% =============================== FUNCTIONS ===============================
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

function yi = safe_interp(x, y, xi)
    x = x(:); y = y(:); xi = xi(:);
    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    if numel(x) < 2
        yi = nan(size(xi));
        return;
    end
    [x, idx] = sort(x); y = y(idx);
    yi = interp1(x, y, xi, 'linear', 'extrap');
end

function [Pxx, N, cfg_out] = welch_psd_spatial(x, fs_spatial, nfft_user, win_frac, ovlp_frac, remove_dc)
    x = x(:);
    if remove_dc
        x = x - mean(x,'omitnan');
    end
    Npts = numel(x);
    if isempty(nfft_user)
        nfft = max(1024, 2^nextpow2(min(8192, Npts)));
    else
        nfft = max(256, nfft_user);
    end
    wlen  = max(64, round(win_frac * nfft));
    wlen  = min(wlen, Npts);
    if mod(wlen,2)==0, wlen = wlen-1; end
    nover = max(0, round(ovlp_frac * wlen));
    win   = hamming(wlen);
    [Pxx, N] = pwelch(x, win, nover, nfft, fs_spatial, 'onesided');
    cfg_out = struct('nfft',nfft,'wlen',wlen,'nover',nover,'remove_dc',logical(remove_dc));
end

function y2 = second_derivative(y, dx)
    y = y(:);
    if numel(y) < 5
        y2 = zeros(size(y));
        return;
    end
    % Central-difference second derivative
    y2 = zeros(size(y));
    for i=2:numel(y)-1
        y2(i) = (y(i+1) - 2*y(i) + y(i-1)) / (dx^2);
    end
    y2(1) = y2(2);
    y2(end) = y2(end-1);
end

function idx = pick_zoom_centers(x, strength, K, min_sep_m)
    % Pick up to K indices where "strength" is large, separated by min_sep_m.
    x = x(:); strength = strength(:);
    if isempty(strength) || all(~isfinite(strength)), idx = []; return; end

    % Normalize and rank candidates
    s = strength; s(~isfinite(s)) = 0;
    [~, order] = sort(s, 'descend');

    % Greedy selection with separation constraint
    idx = [];
    for ii = 1:numel(order)
        if numel(idx) >= K, break; end
        i = order(ii);
        if isempty(idx)
            idx = i;
        else
            if all(abs(x(i) - x(idx)) >= min_sep_m)
                idx(end+1) = i; %#ok<AGROW>
            end
        end
    end

    % Sort by position
    idx = sort(idx);
end

function [a1, a0, r] = simple_linreg(x, y)
    x = x(:); y = y(:);
    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    if numel(x) < 3
        a1 = NaN; a0 = NaN; r = NaN; return;
    end
    X = [x ones(size(x))];
    beta = X \ y;
    a1 = beta(1); a0 = beta(2);
    C = corrcoef(x, y);
    if all(size(C)>=2), r = C(1,2); else r = NaN; end
end

function add_vline(xpos, style, labelStr)
    if isempty(xpos) || ~isfinite(xpos), return; end
    yl = ylim; hold on;
    plot([xpos xpos], yl, style, 'Color',[0.2 0.2 0.2], 'LineWidth', 1.0);
    text(xpos, yl(2), sprintf('  %s (%.3f)', labelStr, xpos), ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'Rotation',90,'FontSize',9,'Color',[0.2 0.2 0.2]);
end