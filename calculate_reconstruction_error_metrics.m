
%% 7) Calculate Quantitative Error Metrics vs. Ground Truth
% -------------------------------------------------------------------------
% PURPOSE:
% This script performs a quantitative comparison between the final
% reconstructed road profile and the original ground truth profile. It
% calculates a suite of standard error metrics to objectively measure the
% accuracy of the reconstruction.
%
% CORE LOGIC (How it works):
% 1.  Loads the reconstructed spatial profile ('recon_spatial.mat') and the
%     configuration file ('recon_cfg.mat'), which contains the ground truth
%     road profile.
%
% 2.  Places both the reconstructed and true profiles onto a common spatial
%     grid by interpolating one onto the other. This ensures a direct,
%     point-for-point comparison is possible.
%
% 3.  Calculates pointwise error metrics by comparing the two profiles at
%     each location. This includes:
%     - Bias (the average error)
%     - RMSE (Root Mean Square Error)
%     - MAE (Mean Absolute Error)
%     - Pearson Correlation Coefficient (r)
%
% 4.  Calculates spectral error metrics by comparing the Power Spectral
%     Density (PSD) of the two profiles. This measures how well the
%     reconstruction captures the road's roughness at different spatial
%     frequencies (wavelengths). This includes the overall Relative
%     Spectral Error (RSE) and errors within specific frequency bands.
%
% 5.  Saves all calculated metrics into a .mat file and a summary .txt file.
%
% INPUTS:
% - '02_Reconstruction Inputs/recon_cfg.mat': For the ground truth profile.
% - '04_Spatial Mapping/recon_spatial.mat': For the reconstructed profile.
%
% OUTPUTS:
% - '07_Comparison Metrics/recon_compare_metrics.mat': Contains all calculated
%   error metrics.
% - '07_Comparison Metrics/recon_compare_metrics.txt': A human-readable
%   summary of the key metrics.
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
save_folder = fullfile('00_Outputs', '04_Validation', 'ErrorMetrics');
fig_dir     = fullfile(save_folder,'figs');

% (3) Which reconstructed series to compare?
%     'filtered' -> use zr_x_filt (recommended; paper low-pass policy)
%     'raw'      -> use zr_x_raw  (pre-filter)
recon_series = 'filtered';  % 'filtered' | 'raw'

% (4) Spatial-frequency bands for spectral error (cycles/m)
%     Defaults are broad & practical; they will be clipped to Nyquist automatically.
bands_cyc_per_m = [ 0.00  0.30;
                    0.30  0.80;
                    0.80  inf ];

% (5) Spectral analysis (Welch PSD) settings (spatial domain)
welch_nfft_points  = [];     % [] => auto from length; else integer (power of 2 is fine)
welch_window_frac  = 0.50;   % window length as fraction of N (0.3–0.6 typical)
welch_overlap_frac = 0.50;   % 50% overlap

% (6) Metrics options
remove_dc_for_rmse = false;  % if true, subtract mean before RMSE/MAE; keep false to reflect bias separately
remove_dc_for_psd  = true;   % recommended: spectra should ignore DC
use_common_support = true;   % use overlapping x-range only (recommended)

% (7) Plots & export
make_plots = true;
export_png = true;
% =====================================================================

%% ------------------------- Locate & load inputs -------------------------
cfg_path = '';
for k = 1:numel(cfg_candidates)
    if exist(cfg_candidates{k},'file') == 2
        cfg_path = cfg_candidates{k};
        break;
    end
end
assert(~isempty(cfg_path), 'Missing recon_cfg.mat in expected locations.');

sp_path = '';
for k = 1:numel(sp_candidates)
    if exist(sp_candidates{k},'file') == 2, sp_path = sp_candidates{k}; break; end
end
assert(~isempty(sp_path), 'Missing recon_spatial.mat in expected locations.');

% Make output folders
if ~exist(save_folder,'dir'), mkdir(save_folder); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

C  = load(cfg_path); cfg = C.cfg;
SP = load(sp_path);  sp  = SP.sp;

% Ground truth (required)
assert(isfield(cfg,'ground') && isfield(cfg.ground,'x') && ~isempty(cfg.ground.x), ...
    'Ground-truth profile not found in cfg.ground.x. Provide ground truth to compare.');
x_true  = cfg.ground.x(:);
zr_true = cfg.ground.zr(:);

% Reconstructed series
x_rec   = sp.x(:);
switch lower(recon_series)
    case 'filtered', zr_rec = sp.zr_x_filt(:);
    case 'raw',      zr_rec = sp.zr_x_raw(:);
    otherwise, error('Unknown recon_series="%s"', recon_series);
end

% Restrict to common support (optional)
if use_common_support
    xmin = max(min(x_true), min(x_rec));
    xmax = min(max(x_true), max(x_rec));
    mT = (x_true>=xmin) & (x_true<=xmax);
    mR = (x_rec >=xmin) & (x_rec <=xmax);
    xT = x_true(mT); yT = zr_true(mT);
    xR = x_rec(mR);  yR = zr_rec(mR);
else
    xT = x_true; yT = zr_true;
    xR = x_rec;  yR = zr_rec;
end

% Interpolate reconstruction to truth grid (reference = ground truth)
zr_rec_on_true = safe_interp(xR, yR, xT);

% Spatial sampling info from reconstruction (assume similar to truth)
dx         = mean(diff(xT),'omitnan');
fs_spatial = 1 / max(dx, eps);     % samples per meter
nyq_sp     = fs_spatial / 2;       % cycles per meter

% Clip band edges to [0, nyq]
bands = bands_cyc_per_m;
bands(:,1) = max(bands(:,1), 0);
bands(:,2) = min(bands(:,2), nyq_sp);
bands(bands(:,2) <= bands(:,1), :) = [];  % remove empty bands

%% ------------------------- Pointwise metrics -------------------------
y_true = yT(:);
y_rec  = zr_rec_on_true(:);

% Bias & DC handling
bias_val  = mean(y_rec - y_true, 'omitnan');
if remove_dc_for_rmse
    e = (y_rec - mean(y_rec,'omitnan')) - (y_true - mean(y_true,'omitnan'));
else
    e = y_rec - y_true;
end

rmse    = sqrt(mean(e.^2, 'omitnan'));
mae     = mean(abs(e), 'omitnan');
r_val   = robust_corr(y_true, y_rec);
nrmse_r = rmse / (max(y_true)-min(y_true) + eps);           % vs range
nrmse_s = rmse / (std(y_true,0,'omitnan') + eps);           % vs std

% Endpoint drift / net elevation change
d_true = y_true(end) - y_true(1);
d_rec  = y_rec(end)  - y_rec(1);
endpoint_delta      = d_rec - d_true;
endpoint_delta_norm = endpoint_delta / (std(y_true,0,'omitnan') + eps);

%% ------------------------- Spectral metrics -------------------------
% Welch PSDs (spatial)
[Ptrue, Ntrue, welch_cfg] = welch_psd_spatial(y_true, fs_spatial, ...
                        welch_nfft_points, welch_window_frac, welch_overlap_frac, remove_dc_for_psd);
[Prec,  Nrec]             = welch_psd_spatial(y_rec, fs_spatial, ...
                        welch_nfft_points, welch_window_frac, welch_overlap_frac, remove_dc_for_psd);

% Interpolate PSDs to a common spatial frequency grid
n_common = linspace(0, min([Ntrue(end), Nrec(end), nyq_sp]), max(numel(Ntrue), numel(Nrec)))';
Ptrue_c  = safe_interp(Ntrue, Ptrue, n_common);
Prec_c   = safe_interp(Nrec,  Prec,  n_common);

% Overall relative spectral error (RSE)
rse_all = sqrt( sum((Prec_c - Ptrue_c).^2) / (sum(Ptrue_c.^2) + eps) );

% Bandwise integrated power errors and RSE
bands_metrics = struct('n1',{},'n2',{},'power_true',{},'power_rec',{},'power_ratio',{},'power_rel_err',{},'rse',{});
for k = 1:size(bands,1)
    n1 = bands(k,1);
    n2 = bands(k,2);
    m  = (n_common>=n1) & (n_common<n2);
    Pt = Ptrue_c(m);
    Pr = Prec_c(m);
    pow_t = trapz(n_common(m), Pt);
    pow_r = trapz(n_common(m), Pr);
    ratio = pow_r / (pow_t + eps);
    relerr= abs(pow_r - pow_t) / (pow_t + eps);
    rse_k = sqrt( sum((Pr - Pt).^2) / (sum(Pt.^2) + eps) );

    bands_metrics(k).n1 = n1;
    bands_metrics(k).n2 = n2;
    bands_metrics(k).power_true = pow_t;
    bands_metrics(k).power_rec  = pow_r;
    bands_metrics(k).power_ratio= ratio;
    bands_metrics(k).power_rel_err = relerr;
    bands_metrics(k).rse = rse_k;
end

%% ------------------------- Save metrics -------------------------
M = struct();
M.info.recon_series         = recon_series;
M.info.use_common_support   = use_common_support;
M.info.fs_spatial           = fs_spatial;
M.info.nyquist_cyc_per_m    = nyq_sp;
M.info.bands_cyc_per_m      = bands;

M.pointwise.bias_m          = bias_val;
M.pointwise.rmse_m          = rmse;
M.pointwise.mae_m           = mae;
M.pointwise.pearson_r       = r_val;
M.pointwise.nrmse_range     = nrmse_r;
M.pointwise.nrmse_std       = nrmse_s;
M.pointwise.endpoint_delta_m = endpoint_delta;
M.pointwise.endpoint_delta_norm = endpoint_delta_norm;

M.spectral.rse_all          = rse_all;
M.spectral.welch_cfg        = welch_cfg;
M.spectral.n_common         = n_common;
M.spectral.P_true_c         = Ptrue_c;
M.spectral.P_rec_c          = Prec_c;
M.spectral.bands            = bands_metrics;

save(fullfile(save_folder,'recon_compare_metrics.mat'), 'M','-v7.3');
fprintf('Saved: %s\n', fullfile(save_folder,'recon_compare_metrics.mat'));

% Text summary
txt_path = fullfile(save_folder,'recon_compare_metrics.txt');
fid = fopen(txt_path,'w');
if fid>0
    fprintf(fid, 'Reconstruction vs Truth — Metrics Summary\n');
    fprintf(fid, '========================================\n');
    fprintf(fid, 'Series compared : %s\n', recon_series);
    fprintf(fid, 'Common support  : %d\n', use_common_support);
    fprintf(fid, 'Spatial fs (1/m): %.5f   | Nyquist: %.5f cyc/m\n', fs_spatial, nyq_sp);
    fprintf(fid, '\n--- Pointwise ---\n');
    fprintf(fid, 'Bias (m)        : %+ .6e\n', M.pointwise.bias_m);
    fprintf(fid, 'RMSE (m)        : %.6e\n',  M.pointwise.rmse_m);
    fprintf(fid, 'MAE  (m)        : %.6e\n',  M.pointwise.mae_m);
    fprintf(fid, 'Pearson r       : %.5f\n',   M.pointwise.pearson_r);
    fprintf(fid, 'NRMSE (range)   : %.5f\n',   M.pointwise.nrmse_range);
    fprintf(fid, 'NRMSE (std)     : %.5f\n',   M.pointwise.nrmse_std);
    fprintf(fid, 'Endpoint Δ (m)  : %+ .6e\n', M.pointwise.endpoint_delta_m);
    fprintf(fid, 'Endpoint Δ /σ   : %.5f\n',   M.pointwise.endpoint_delta_norm);
    fprintf(fid, '\n--- Spectral ---\n');
    fprintf(fid, 'RSE (overall)   : %.5f\n',   M.spectral.rse_all);
    for k=1:numel(bands_metrics)
        fprintf(fid, 'Band [%.3f, %.3f) cyc/m: Power ratio=%.3f, RelErr=%.3f, RSE=%.5f\n', ...
            bands_metrics(k).n1, bands_metrics(k).n2, ...
            bands_metrics(k).power_ratio, bands_metrics(k).power_rel_err, bands_metrics(k).rse);
    end
    fclose(fid);
end
fprintf('Saved: %s\n', txt_path);

%% ------------------------- Figures (improved) -------------------------
if make_plots
    % ------- Style / palette -------
    baseFont = 12; axesLW = 0.9; lineLW = 1.8; useMinorGrid = true;
    C.black = [0.10 0.10 0.10];
    C.blue  = [0.00 0.45 0.74];
    C.red   = [0.85 0.33 0.10];
    C.gray  = [0.60 0.60 0.60];
    C.green = [0.20 0.62 0.20];
    C.purp  = [0.49 0.18 0.56];

    set(0,'defaultAxesFontName','Helvetica','defaultAxesFontSize',baseFont);
    set(0,'defaultLineLineWidth',lineLW);
    nyq_local = nyq_sp;

    % ====== Fig 1 — metrics summary panel ======
    f1 = figure('Name','Metrics summary','Color','w','Units','pixels'); f1.Position(3:4)=[860 440];
    axes('Position',[0 0 1 1]); axis off;
    lines = {
        sprintf('Series: %s    | common support: %d', recon_series, use_common_support)
        sprintf('fs_{spatial} = %.4f  (Nyquist = %.4f cyc/m)', fs_spatial, nyq_sp)
        ' '
        sprintf('Bias (m)         : %+ .6e', bias_val)
        sprintf('RMSE (m)         : %.6e',  rmse)
        sprintf('MAE  (m)         : %.6e',  mae)
        sprintf('Pearson r        : %.5f',   r_val)
        sprintf('NRMSE (range)    : %.5f',   nrmse_r)
        sprintf('NRMSE (std)      : %.5f',   nrmse_s)
        sprintf('Endpoint Δ (m)   : %+ .6e', endpoint_delta)
        sprintf('Endpoint Δ /σ    : %.5f',   endpoint_delta_norm)
        ' '
        sprintf('RSE (overall)    : %.5f',   rse_all)
    };
    text(0.06,0.90, lines, 'Units','normalized', 'FontName','Consolas','FontSize',11, 'VerticalAlignment','top');
    if export_png, save_fig_safe(f1, fullfile(fig_dir,'13_metrics_summary.png'), 180); end

    % ====== Fig 2 — overlay + residual ======
    f2 = figure('Name','Aligned overlay + residual','Color','w','Units','pixels'); f2.Position(3:4)=[1000 700];
    tlo2 = tiledlayout(f2,2,1,'Padding','compact','TileSpacing','compact');

    ax21 = nexttile(tlo2,1); hold(ax21,'on'); set(ax21,'LineWidth',axesLW);
    plot(ax21, xT, y_true, '-', 'Color', C.black, 'DisplayName','True');
    plot(ax21, xT, y_rec,  '-', 'Color', C.blue,  'DisplayName','Recon');
    grid(ax21,'on'); if useMinorGrid, grid(ax21,'minor'); end
    ylabel(ax21,'z_r (m)');
    [ymin,ymax] = bounds([y_true; y_rec]);
    pad = 0.05 * max(eps, ymax - ymin);
    ylim(ax21, [ymin-pad, ymax+pad]);
    legend(ax21,'Location','best','Box','off');
    title(ax21, sprintf('Overlay — RMSE=%.4g m, NRMSE(range)=%.2f%%, r=%.3f', rmse, 100*nrmse_r, r_val));

    ax22 = nexttile(tlo2,2); hold(ax22,'on'); set(ax22,'LineWidth',axesLW);
    plot(ax22, xT, e, '-', 'Color', C.red);
    yline(ax22, 0, ':', 'Color', C.gray);
    grid(ax22,'on'); if useMinorGrid, grid(ax22,'minor'); end
    xlabel(ax22,'Distance x (m)'); ylabel(ax22,'Residual (m)');
    ylim(ax22, symmetric_limits(e));
    title(ax22,'Residual (Recon − True) — balanced y-limits');
    linkaxes([ax21, ax22],'x');
    if export_png, save_fig_safe(f2, fullfile(fig_dir,'14_overlay_and_residual.png'), 300); end

    % ====== Fig 3 — scatter (True vs Recon) with fit & 1:1 ======
    f3 = figure('Name','Scatter True vs Recon','Color','w','Units','pixels'); f3.Position(3:4)=[720 560];
    ax3 = axes(f3); hold(ax3,'on'); set(ax3,'LineWidth',axesLW);
    m = isfinite(y_true) & isfinite(y_rec);
    scatter(ax3, y_true(m), y_rec(m), 8, 'MarkerEdgeColor', C.blue, 'MarkerFaceColor', C.blue,...
        'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.25);
    % 1:1
    yy = [min(y_true(m)); max(y_true(m))]; plot(ax3, yy, yy, '-', 'Color', C.gray, 'DisplayName','1:1');
    % Fit
    if nnz(m) > 2
        p = polyfit(y_true(m), y_rec(m), 1); yfit = polyval(p, yy);
        plot(ax3, yy, yfit, '--', 'Color', C.purp, 'DisplayName',sprintf('fit: y=%.3fx%+ .3g', p(1), p(2)));
        r2 = max(0, corr(y_true(m), y_rec(m),'Rows','pairwise')^2);
    else
        p = [NaN NaN]; r2 = NaN;
    end
    grid(ax3,'on'); if useMinorGrid, grid(ax3,'minor'); end
    xlabel(ax3,'True z_r (m)'); ylabel(ax3,'Recon z_r (m)');
    title(ax3, sprintf('True vs Recon — slope=%.3f, intercept=%+.3g,  R^2=%.3f', p(1), p(2), r2));
    legend(ax3,'Location','best');
    if export_png, save_fig_safe(f3, fullfile(fig_dir,'15_scatter_true_vs_recon.png'), 300); end

    % ====== Fig 4 — PSD (common grid) with band shading ======
    % Use your already-interpolated common grid
    kpos = find(n_common>0,1,'first'); fplt = n_common(kpos:end);
    Ptrue_db = 10*log10(Ptrue_c(kpos:end) + eps);
    Prec_db  = 10*log10(Prec_c(kpos:end)  + eps);

    f4 = figure('Name','PSD compare','Color','w','Units','pixels'); f4.Position(3:4)=[1000 600];
    tlo4 = tiledlayout(f4,1,1,'Padding','compact','TileSpacing','compact');
    ax4 = nexttile(tlo4); hold(ax4,'on'); set(ax4,'XScale','log','LineWidth',axesLW);
    plot(ax4, fplt, Ptrue_db, '-', 'Color', C.black, 'DisplayName','True PSD');
    plot(ax4, fplt, Prec_db,  '-', 'Color', C.blue,  'DisplayName','Recon PSD');
    grid(ax4,'on'); if useMinorGrid, grid(ax4,'minor'); end
    xlabel(ax4,'Spatial frequency (cycles/m)'); ylabel(ax4,'PSD (dB re m^2/(cycles/m))');
    title(ax4,'Spatial PSD — True vs Recon');
    drawnow; yl = ylim(ax4);
    shade_bands_psd(ax4, fplt, yl, bands, nyq_local);
    legend(ax4,'Location','best','Box','off');
    if export_png, save_fig_safe(f4, fullfile(fig_dir,'16_psd_compare.png'), 300); end

    % ====== Fig 5 — PSD ratio (Recon/True) in dB ======
    ratio_db = 10*log10((Prec_c(kpos:end)+eps) ./ (Ptrue_c(kpos:end)+eps));
    f5 = figure('Name','PSD ratio','Color','w','Units','pixels'); f5.Position(3:4)=[1000 420];
    ax5 = axes(f5); hold(ax5,'on'); set(ax5,'XScale','log','LineWidth',axesLW);
    plot(ax5, fplt, ratio_db, 'Color', C.purp);
    yline(ax5, 0, '-', 'Color', C.gray, 'DisplayName','0 dB');
    yline(ax5, 3, ':', 'Color', C.gray, 'HandleVisibility','off');
    yline(ax5,-3, ':', 'Color', C.gray, 'HandleVisibility','off');
    grid(ax5,'on'); if useMinorGrid, grid(ax5,'minor'); end
    xlabel(ax5,'Spatial frequency (cycles/m)'); ylabel(ax5,'Recon/True PSD (dB)');
    title(ax5,'PSD mismatch (positive = overestimate)');
    drawnow; yl = ylim(ax5);
    shade_bands_psd(ax5, fplt, yl, bands, nyq_local, 0.05);
    if export_png, save_fig_safe(f5, fullfile(fig_dir,'17_psd_ratio.png'), 300); end

    % ====== Fig 6 — Band metrics bars ======
    if ~isempty(bands_metrics)
        br = arrayfun(@(b) b.power_ratio,   bands_metrics);
        be = arrayfun(@(b) b.power_rel_err, bands_metrics);
        rr = arrayfun(@(b) b.rse,           bands_metrics);
        labels = arrayfun(@(b) sprintf('[%.2f,%.2f)', b.n1, b.n2), bands_metrics, 'UniformOutput', false);

        f6 = figure('Name','Spectral band metrics','Color','w','Units','pixels'); f6.Position(3:4)=[1000 640];
        tlo6 = tiledlayout(f6,3,1,'Padding','compact','TileSpacing','compact');

        ax61 = nexttile(tlo6,1); hold(ax61,'on'); set(ax61,'LineWidth',axesLW);
        bh1 = bar(ax61, br, 'FaceColor', C.blue); grid(ax61,'on'); if useMinorGrid, grid(ax61,'minor'); end
        ylabel(ax61,'Power ratio'); title(ax61,'rec/true'); xticks(ax61,1:numel(br)); xticklabels(ax61,labels);
        label_bars(ax61, br);

        ax62 = nexttile(tlo6,2); hold(ax62,'on'); set(ax62,'LineWidth',axesLW);
        bh2 = bar(ax62, be, 'FaceColor', C.purp); grid(ax62,'on'); if useMinorGrid, grid(ax62,'minor'); end
        ylabel(ax62,'|ΔP| / P_{true}'); title(ax62,'Relative power error'); xticks(ax62,1:numel(be)); xticklabels(ax62,labels);
        label_bars(ax62, be);

        ax63 = nexttile(tlo6,3); hold(ax63,'on'); set(ax63,'LineWidth',axesLW);
        bh3 = bar(ax63, rr, 'FaceColor', C.red); grid(ax63,'on'); if useMinorGrid, grid(ax63,'minor'); end
        ylabel(ax63,'RSE'); xlabel(ax63,'Bands (cycles/m)'); title(ax63,'Relative spectral error');
        xticks(ax63,1:numel(rr)); xticklabels(ax63,labels);
        label_bars(ax63, rr);

        if export_png, save_fig_safe(f6, fullfile(fig_dir,'18_spectral_error_bands.png'), 300); end
    end

    % ====== Fig 7 — Residual histogram ======
    f7 = figure('Name','Residual histogram','Color','w','Units','pixels'); f7.Position(3:4)=[760 480];
    ax7 = axes(f7); hold(ax7,'on'); set(ax7,'LineWidth',axesLW);
    histogram(ax7, e, 50, 'Normalization','pdf', 'FaceColor', C.red, 'FaceAlpha',0.6, 'EdgeColor','none');
    mu = mean(e,'omitnan'); sg = std(e,'omitnan');
    xline(ax7, 0,  ':', 'Color', C.gray);
    xline(ax7, mu, '-', 'Color', C.black, 'DisplayName','mean');
    xline(ax7, mu+sg,  ':', 'Color', C.gray, 'HandleVisibility','off');
    xline(ax7, mu-sg,  ':', 'Color', C.gray, 'HandleVisibility','off');
    xline(ax7, mu+2*sg,':', 'Color', C.gray, 'HandleVisibility','off');
    xline(ax7, mu-2*sg,':', 'Color', C.gray, 'HandleVisibility','off');
    grid(ax7,'on'); if useMinorGrid, grid(ax7,'minor'); end
    xlabel(ax7,'Residual (m)'); ylabel(ax7,'PDF');
    title(ax7, sprintf('Residual distribution — \\mu=%.4g m, \\sigma=%.4g m', mu, sg));
    if export_png, save_fig_safe(f7, fullfile(fig_dir,'19_residual_hist.png'), 300); end
end

%% =============================== FUNCTIONS ===============================
function [Pxx, N, cfg_out] = welch_psd_spatial(x, fs_spatial, nfft_user, win_frac, ovlp_frac, remove_dc)
    % Welch PSD in spatial domain (x: profile vs distance).
    % Returns Pxx (power spectral density), N (cycles/m), and config used.
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
    if mod(wlen,2)==0, wlen = wlen-1; end  % odd window length is fine
    nover = max(0, round(ovlp_frac * wlen));
    win   = hamming(wlen);
    [Pxx, N] = pwelch(x, win, nover, nfft, fs_spatial, 'onesided');

    cfg_out = struct('nfft',nfft,'wlen',wlen,'nover',nover,'remove_dc',logical(remove_dc));
end

function r = robust_corr(a, b)
    a = a(:); b = b(:);
    m = isfinite(a) & isfinite(b);
    a = a(m); b = b(m);
    if numel(a) < 3
        r = NaN; return;
    end
    C = corrcoef(a, b);
    if all(size(C) >= 2) && isfinite(C(1,2))
        r = C(1,2);
    else
        r = NaN;
    end
end

function yi = safe_interp(x, y, xi)
    x = x(:); y = y(:); xi = xi(:);
    % Remove NaNs in source
    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    if numel(x) < 2
        yi = nan(size(xi));
        return;
    end
    % Ensure monotonic for interp1
    [x, idx] = sort(x); y = y(idx);
    yi = interp1(x, y, xi, 'linear', 'extrap');
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

function save_fig_safe(fig, filepath, dpi)
% Export without transparency requests; fallback to print if needed.
    try
        exportgraphics(fig, filepath, 'Resolution', dpi);
    catch
        [p,f,e] = fileparts(filepath);
        if isempty(e), e = '.png'; end
        switch lower(e)
            case '.png', print(fig, fullfile(p,[f e]), '-dpng',  sprintf('-r%d',dpi));
            case '.tif', print(fig, fullfile(p,[f e]), '-dtiff', sprintf('-r%d',dpi));
            case '.pdf'
                set(fig,'PaperPositionMode','auto'); 
                print(fig, fullfile(p,[f e]), '-dpdf','-vector');
            otherwise,  print(fig, fullfile(p,[f '.png']), '-dpng', sprintf('-r%d',dpi));
        end
    end
    try savefig(fig, [filepath(1:end-4) '.fig']); catch, end
end

function shade_bands_psd(ax, f, ylims, bands, nyq, alpha)
% Shade low/mid/high bands for PSD-like plots on log-x.
    if nargin<6 || isempty(alpha), alpha = 0.08; end
    hold(ax,'on');
    Cb = [0.80 0.90 1.00];  % low
    Cm = [0.92 0.85 0.96];  % mid
    Ch = [1.00 0.88 0.88];  % high
    xmin = f(1); xmax = f(end);
    for k=1:size(bands,1)
        x1 = max(xmin, max(0, bands(k,1)));
        x2 = min(xmax, min(nyq*0.999, bands(k,2)));
        if x2 > x1
            col = Cb; if k==2, col=Cm; elseif k==3, col=Ch; end
            patch('Parent',ax,'XData',[x1 x2 x2 x1],'YData',[ylims(1) ylims(1) ylims(2) ylims(2)], ...
                  'FaceColor',col,'FaceAlpha',alpha,'EdgeColor','none','HandleVisibility','off');
        end
    end
    % Label separators
    for k=1:size(bands,1)-1
        xline(ax, bands(k,2), ':', 'Color',[0.4 0.4 0.4], 'HandleVisibility','off');
    end
end

function yL = symmetric_limits(y)
% Symmetric balanced limits about zero using robust max (98th percentile)
    y = y(isfinite(y));
    if isempty(y), yL = [-1 1]; return; end
    m = prctile(abs(y),98);
    m = max(m, max(abs(y))*0.5); % avoid too tight if heavy tails
    yL = 1.05*[-m m];
end

function label_bars(ax, vals)
% Add value labels above bars
    for i=1:numel(vals)
        text(ax, i, vals(i), sprintf('%.3g', vals(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end