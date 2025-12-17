
%% 4) Profile Registration and Metrics
% -----------------------------------------------------------------------------
% PURPOSE
%   Align reconstructed vs. true road profiles in the spatial domain and
%   compute accuracy metrics (RMSE/NRMSE/MAE/Corr/SSIM) + band-limited errors.
%   This version AUTO-DETECTS the reconstructed profile regardless of the
%   chosen reconstruction path:
%     1) 04_Spatial Mapping and LP/recon_spatial.mat   (struct 'sp')  [PRIMARY]
%     2) 03B_TF Inversion/recon_spatial.mat            (struct 'sp')  [FALLBACK]
%     3) 03A_Two-Step ODE Inversion/recon_time_domain.mat + cfg       [BUILD]
%   It also supports legacy 'pkg.spatial.*' if encountered.
%
% INPUTS (auto-discovered):
%   • Reconstructed spatial profile (preferred):
%       sp.x [m], sp.zr_x_filt or sp.zr_x_raw [m], sp.fs_spatial [samples/m]
%   • OR legacy:
%       pkg.spatial.x, pkg.spatial.zr_filt/zr/zr_raw, pkg.spatial.fs_spatial
%   • True profile (sim):
%       ./01_Quarter Car Modeling/simulation_results.mat
%       (accepts {x_true,zr_true} or {t,v,zr_true}, resampled to recon grid)
%
% OUTPUTS:
%   ./07_Registration_Metrics/
%     aligned_pair.mat            (struct 'align' with x, true, recon, resid, lag)
%     registration_metrics.mat    (struct 'metrics' + 'bands')
%     registration_metrics.csv    (flat table)
%     figs/*                      (informative figures)
%
% NOTES:
%   • Plots and filters are Nyquist-safe.
%   • Spatial bands are given in cycles/m (adjust below as needed).
% -----------------------------------------------------------------------------

%% =========================
%  [1] HYPERPARAMETERS
%  =========================
clear; clc; close all;

% Output dirs
out_dir  = fullfile('00_Outputs', '04_Validation', 'ErrorMetrics');
fig_dir  = fullfile(out_dir,'figs');

if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Where to look for reconstructed results (in order)
spatial_candidates = { ...
    fullfile('00_Outputs', '03_ReconstructionCore', 'SpatialMapping', 'recon_spatial.mat'), ...        % From Path A
    fullfile('00_Outputs', '03_ReconstructionCore', 'PathB_FrequencyDomain', 'recon_spatial.mat'), ... % From Path B
};

% Legacy 'pkg' fallbacks (if you used that schema before)
pkg_candidates = { ...
    fullfile('00_Outputs', '06_FinalExports', 'ConsolidatedPackage', 'reconstructed_profile.mat'), ...
};

% If only time-domain exists (Path A), we can build a spatial profile on the fly:
pathA_time_candidate = fullfile('00_Outputs', '03_ReconstructionCore', 'PathA_TimeDomain', 'recon_time_domain.mat');
cfg_candidates       = { ...
    fullfile('00_Outputs', '02_ReconstructionSetup', 'Config', 'recon_cfg.mat') ...
};

% Simulation truth
sim_mat_path = fullfile('00_Outputs', '01_Simulation', 'SimulationData', 'simulation_results.mat');

% Spatial bands (cycles/m) for error decomposition
bands.c_low.max  = 0.05;     % long waves
bands.c_mid.min  = 0.05;
bands.c_mid.max  = 0.30;
bands.c_high.min = 0.30;     % short waves

% Welch PSD
psd_nfft_max = 8192;

% XCorr alignment search (fraction of length)
max_shift_fraction = 0.10;

% Plot export
dpi = 300;

% ---- Display controls ----
plot_overlay_metrics   = false;   % draw metrics box on the figure (OFF)
print_metrics_to_cmd   = true;    % print metrics to Command Window (ON)

%% =========================
%  [2] LOAD RECONSTRUCTED SPATIAL PROFILE (sp or pkg)
%  =========================
recon_source = '';
xq = [];
zr_recon = [];
fs_spatial = [];
dx = [];

% Try 'sp' first
sp_path = first_existing(spatial_candidates);
if ~isempty(sp_path)
    S = load(sp_path);
    assert(isfield(S,'sp'), 'Expected struct ''sp'' in %s', sp_path);
    sp = S.sp;

    % Prefer filtered profile if available
    zr_recon = choose_first_field(sp, {'zr_x_filt','zr_x','zr_x_raw'});
    xq        = sp.x(:);
    fs_spatial= get_first(sp, {'fs_spatial'}, []);
    if isempty(fs_spatial)
        dx = mean(diff(xq),'omitnan'); fs_spatial = 1/max(dx,eps);
    else
        dx = 1/fs_spatial;
    end
    recon_source = sp_path;

else
    % Try legacy 'pkg.spatial'
    pkg_path = first_existing(pkg_candidates);
    if ~isempty(pkg_path)
        P = load(pkg_path);
        pkg = [];
        if isfield(P,'pkg'), pkg = P.pkg; end
        if isempty(pkg) && isfield(P,'recon') && isfield(P.recon,'pkg')
            pkg = P.recon.pkg;
        end
        assert(~isempty(pkg) && isfield(pkg,'spatial'), 'No pkg.spatial found in %s', pkg_path);

        zr_recon    = choose_first_field(pkg.spatial, {'zr_filt','zr','zr_raw'});
        xq          = pkg.spatial.x(:);
        fs_spatial  = get_first(pkg.spatial, {'fs_spatial'}, []);
        if isempty(fs_spatial)
            dx = mean(diff(xq),'omitnan'); fs_spatial = 1/max(dx,eps);
        else
            dx = 1/fs_spatial;
        end
        recon_source = pkg_path;

    else
        % Last chance: build spatial from Path A time-domain + cfg
        assert(exist(pathA_time_candidate,'file')==2, ...
            ['Reconstructed spatial profile not found. Neither ''sp'' nor ''pkg'' found.\n' ...
             'Also missing shared spatial file. Run step 4 mapping OR Path B publish.\n' ...
             'Expected at least: %s'], pathA_time_candidate);

        RT  = load(pathA_time_candidate);  assert(isfield(RT,'recon'),'Missing recon.* in time-domain mat');
        cfgp= first_existing(cfg_candidates);  assert(~isempty(cfgp),'Could not locate recon_cfg.mat');
        C   = load(cfgp);  assert(isfield(C,'cfg'),'recon_cfg.mat must contain cfg');

        t     = RT.recon.t(:);
        zr_t  = RT.recon.zr_t(:);
        V     = C.cfg.V;
        x_t   = get_first(C.cfg, {'x_t'}, []);
        if isempty(x_t), x_t = (t - t(1))*V; end

        % Build a uniform spatial grid over the support of x_t
        xq  = x_t(:);
        zr_recon = interp1(x_t, zr_t, xq, 'linear','extrap');
        zr_recon(~isfinite(zr_recon)) = 0;

        dx = mean(diff(xq),'omitnan');
        fs_spatial = 1/max(dx,eps);
        recon_source = sprintf('%s + %s (built on-the-fly)', pathA_time_candidate, cfgp);
    end
end

% Detrend recon for fair correlation (keep physical mean for plots)
zr_recon = zr_recon(:);
zr_recon_dc = zr_recon - mean(zr_recon,'omitnan');

%% =========================
%  [3] LOAD TRUE PROFILE & RESAMPLE TO xq
%  =========================
assert(exist(sim_mat_path,'file')==2, 'True-profile source not found: %s', sim_mat_path);
T = load(sim_mat_path);

[x_true, zr_true] = extract_truth(T);
[x_true, uidx] = unique(x_true,'stable'); zr_true = zr_true(uidx);
zr_true_rs = interp1(x_true, zr_true, xq, 'linear','extrap');
zr_true_rs = zr_true_rs(:);
zr_true_dc = zr_true_rs - mean(zr_true_rs,'omitnan');

%% =========================
%  [4] ALIGNMENT VIA CROSS-CORR (SPATIAL)
%  =========================
N = numel(xq);
maxK = max(1, round(max_shift_fraction * N));

[r, lags] = xcorr(zr_true_dc, zr_recon_dc, maxK, 'coeff');
[peakR, idxMax] = max(r);
bestLag = lags(idxMax);
lag_m   = bestLag / fs_spatial;

% Apply shift (keep only overlap)
if bestLag >= 0
    recon_al = zr_recon(1+bestLag:end);
    true_al  = zr_true_rs(1:end-bestLag);
    x_al     = xq(1:end-bestLag);
else
    k = -bestLag;
    recon_al = zr_recon(1:end-k);
    true_al  = zr_true_rs(1+k:end);
    x_al     = xq(1:end-k);
end

resid = true_al - recon_al;

%% =========================
%  [5] METRICS
%  =========================
rmse  = sqrt(mean(resid.^2));
mae   = mean(abs(resid));
rngT  = range(true_al);
stdT  = std(true_al);
nrmse_range = rmse / max(rngT, eps);
nrmse_std   = rmse / max(stdT, eps);
Cmat = corrcoef(true_al, recon_al);
cc   = (numel(Cmat)>=4) * Cmat(1,2); if ~cc, cc = NaN; end

ssim1d_val = ssim_1d(recon_al, true_al);

metrics = struct();
metrics.rmse        = rmse;
metrics.mae         = mae;
metrics.nrmse_range = nrmse_range;
metrics.nrmse_std   = nrmse_std;
metrics.corrcoef    = cc;
metrics.ssim1d      = ssim1d_val;
metrics.peak_xcorr  = peakR;
metrics.bestLag_samp= bestLag;
metrics.bestLag_m   = lag_m;
metrics.dx          = 1/fs_spatial;
metrics.fs_spatial  = fs_spatial;
metrics.N_aligned   = numel(x_al);
metrics.recon_source= recon_source;

%% =========================
%  [6] BAND-LIMITED ERRORS
%  =========================
nyq = fs_spatial/2;
clp = @(c) max(eps, min(c, 0.999*nyq));

Wnl = clp(bands.c_low.max) / nyq;
Wnm = [clp(bands.c_mid.min) clp(bands.c_mid.max)] / nyq;
if Wnm(2) <= Wnm(1), Wnm(2) = min(0.999, Wnm(1)+0.05); end
Wnh = clp(bands.c_high.min) / nyq;

[bl, al] = butter(4, Wnl, 'low');
[bm, am] = butter(4, Wnm, 'bandpass');
[bh, ah] = butter(4, Wnh, 'high');

res_L = filtfilt(bl, al, resid);
res_M = filtfilt(bm, am, resid);
res_H = filtfilt(bh, ah, resid);

bands.rmse_low  = sqrt(mean(res_L.^2));
bands.rmse_mid  = sqrt(mean(res_M.^2));
bands.rmse_high = sqrt(mean(res_H.^2));

%% =========================
%  [7] PSD & (optional) coherence
%  =========================
[Ptt, fsp] = welch_psd_safe(true_al,  fs_spatial, psd_nfft_max);
[Prp, ~  ] = welch_psd_safe(recon_al, fs_spatial, psd_nfft_max);

have_mscohere = exist('mscohere','file')==2;
if have_mscohere
    [Cxy, fcoh] = mscohere(true_al, recon_al, [], [], [], fs_spatial);
else
    Cxy = []; fcoh = [];
end

%% =========================
%  [7.5] PRINT METRICS TO COMMAND WINDOW (optional)
%  =========================
if print_metrics_to_cmd
    fprintf('\n[profile_registration_and_metrics]\n');
    fprintf('  Source            : %s\n', recon_source);
    fprintf('  Samples aligned   : %d\n', metrics.N_aligned);
    fprintf('  dx (m)            : %.6f   fs_spatial (samp/m): %.3f\n', metrics.dx, metrics.fs_spatial);
    fprintf('  Best lag          : %.3f m  (%d samples)\n', metrics.bestLag_m, metrics.bestLag_samp);
    fprintf('  RMSE              : %.4f m\n', metrics.rmse);
    fprintf('  MAE               : %.4f m\n', metrics.mae);
    fprintf('  NRMSE (range)     : %.2f %%\n', 100*metrics.nrmse_range);
    fprintf('  NRMSE (std)       : %.2f %%\n', 100*metrics.nrmse_std);
    fprintf('  CorrCoef          : %.3f\n', metrics.corrcoef);
    fprintf('  SSIM1D            : %.3f\n', metrics.ssim1d);
    fprintf('  Band RMSEs (m)    : Low=%.4f  Mid=%.4f  High=%.4f\n\n', ...
            bands.rmse_low, bands.rmse_mid, bands.rmse_high);
end

%% =========================
%  [8] SAVE
%  =========================
align = struct();
align.x                = x_al;
align.zr_true_aligned  = true_al;
align.zr_recon_aligned = recon_al;
align.residual         = resid;
align.bestLag_samples  = bestLag;
align.bestLag_meters   = lag_m;
align.peak_xcorr       = peakR;
align.dx               = 1/fs_spatial;
align.fs_spatial       = fs_spatial;
align.recon_source     = recon_source;

save(fullfile(out_dir,'aligned_pair.mat'), 'align');
save(fullfile(out_dir,'registration_metrics.mat'), 'metrics','bands');

Tmetrics = table(rmse, mae, nrmse_range, nrmse_std, cc, ssim1d_val, peakR, bestLag, lag_m, ...
    'VariableNames', {'RMSE','MAE','NRMSE_range','NRMSE_std','CorrCoef','SSIM1D','PeakXCorr','BestLag_samp','BestLag_m'});
writetable(Tmetrics, fullfile(out_dir,'registration_metrics.csv'));

%% =========================
%  [9] PLOTS
%  =========================

% ---- Style presets ----
baseFont = 12; axesLW = 0.9; lineLW = 1.8; useMinorGrid = true;
C.black = [0.10 0.10 0.10];
C.blue  = [0.00 0.45 0.74];
C.red   = [0.85 0.33 0.10];
C.gray  = [0.60 0.60 0.60];
C.green = [0.20 0.62 0.20];
C.purp  = [0.49 0.18 0.56];
set(0,'defaultAxesFontName','Helvetica','defaultAxesFontSize',baseFont);
set(0,'defaultLineLineWidth',lineLW);
nyq_local = metrics.fs_spatial/2;

% ===================== A) Cross-correlation =====================
figA = figure('Name','A_xcorr','Color','w','Units','pixels'); figA.Position(3:4)=[980 420];
tloA = tiledlayout(figA,1,1,'Padding','compact','TileSpacing','compact');
axA = nexttile(tloA); hold(axA,'on'); set(axA,'LineWidth',axesLW);
xlag = lags / metrics.fs_spatial;
plot(axA, xlag, r, 'Color', C.black, 'DisplayName','r(\tau)');
xline(axA, 0, ':', 'Color', C.gray, 'HandleVisibility','off');
xline(axA, lag_m, '--', 'Color', C.purp, 'HandleVisibility','off');
plot(axA, lag_m, peakR, 'o', 'MarkerSize',6, 'LineWidth',1.2, 'Color', C.purp, 'DisplayName','peak');
grid(axA,'on'); if useMinorGrid, grid(axA,'minor'); end
xlabel(axA,'Lag (m)'); ylabel(axA,'Normalized cross-corr');
title(axA, sprintf('Cross-correlation — peak at %.3f m (r=%.3f)', lag_m, peakR));
xlim(axA, [min(xlag) max(xlag)]);
legend(axA,'Location','best'); 
save_fig(figA, fullfile(fig_dir,'A_xcorr.png'), dpi);

% ===================== B) Aligned profiles + residual =====================
figB = figure('Name','B_aligned_and_residual','Color','w','Units','pixels'); figB.Position(3:4)=[980 700];
tloB = tiledlayout(figB,2,1,'Padding','compact','TileSpacing','compact');

% -- Top: profiles
axB1 = nexttile(tloB,1); hold(axB1,'on'); set(axB1,'LineWidth',axesLW);
plot(axB1, align.x, align.zr_true_aligned, '-', 'Color', C.black, 'DisplayName','True');
plot(axB1, align.x, align.zr_recon_aligned,'-', 'Color', C.blue,  'DisplayName','Recon');
grid(axB1,'on'); if useMinorGrid, grid(axB1,'minor'); end
ylabel(axB1,'Profile z_r (m)');
% Robust y-lims with padding
ymin = min([align.zr_true_aligned; align.zr_recon_aligned]); 
ymax = max([align.zr_true_aligned; align.zr_recon_aligned]);
pad = 0.05 * max(eps, ymax - ymin);
ylim(axB1, [ymin - pad, ymax + pad]);
legend(axB1,'Location','best'); 
title(axB1, sprintf('Aligned profiles — RMSE=%.4g m, NRMSE(range)=%.2f%%, \\rho=%.3f', ...
    metrics.rmse, 100*metrics.nrmse_range, metrics.corrcoef));

% -- Bottom: residual
axB2 = nexttile(tloB,2); hold(axB2,'on'); set(axB2,'LineWidth',axesLW);
plot(axB2, align.x, align.residual, '-', 'Color', C.red);
yline(axB2, 0, ':', 'Color', C.gray);
grid(axB2,'on'); if useMinorGrid, grid(axB2,'minor'); end
xlabel(axB2,'Distance x (m)'); ylabel(axB2,'Residual (m)');
ylim(axB2, symmetric_limits(align.residual));
title(axB2,'Residual (True − Recon) — balanced limits');
linkaxes([axB1, axB2],'x');
save_fig(figB, fullfile(fig_dir,'B_aligned_and_residual.png'), dpi);

% ===================== C) PSD (True vs Recon) with band shading ===========
% Avoid f=0 on log-x
kpos = find(fsp>0,1,'first'); 
fplt = fsp(kpos:end);
Ptt_db = 10*log10(Ptt(kpos:end) + eps);
Prp_db = 10*log10(Prp(kpos:end) + eps);

figC = figure('Name','C_psd_compare','Color','w','Units','pixels'); figC.Position(3:4)=[980 620];
tloC = tiledlayout(figC,1,1,'Padding','compact','TileSpacing','compact');
axC = nexttile(tloC); hold(axC,'on'); set(axC,'LineWidth',axesLW,'XScale','log');
plot(axC, fplt, Ptt_db, 'Color', C.black, 'DisplayName','True');
plot(axC, fplt, Prp_db, 'Color', C.blue,  'DisplayName','Recon');
grid(axC,'on'); if useMinorGrid, grid(axC,'minor'); end
xlabel(axC,'Spatial frequency (cycles/m)'); ylabel(axC,'PSD (dB re m^2/(cycles/m))');
title(axC,'Spatial PSD — True vs Recon');
% Shade bands
drawnow; yl = ylim(axC);
shade_bands(axC, fplt, yl, bands, nyq_local);
legend(axC,'Location','best');
save_fig(figC, fullfile(fig_dir,'C_psd_compare.png'), dpi);

% ===================== D) PSD Ratio (Recon/True) in dB ====================
ratio_db = 10*log10((Prp(kpos:end)+eps)./(Ptt(kpos:end)+eps));
figD = figure('Name','D_psd_ratio','Color','w','Units','pixels'); figD.Position(3:4)=[980 420];
tloD = tiledlayout(figD,1,1,'Padding','compact','TileSpacing','compact');
axD = nexttile(tloD); hold(axD,'on'); set(axD,'LineWidth',axesLW,'XScale','log');
plot(axD, fplt, ratio_db, 'Color', C.purp);
yline(axD, 0, '-', 'Color', C.gray, 'DisplayName','0 dB (perfect)');
yline(axD, 3, ':', 'Color', C.gray, 'HandleVisibility','off');
yline(axD, -3, ':', 'Color', C.gray, 'HandleVisibility','off');
grid(axD,'on'); if useMinorGrid, grid(axD,'minor'); end
xlabel(axD,'Spatial frequency (cycles/m)'); ylabel(axD,'Recon/True PSD (dB)');
title(axD,'PSD mismatch (positive = overestimate)');
drawnow; yl = ylim(axD);
shade_bands(axD, fplt, yl, bands, nyq_local, 0.05); % lighter shade
save_fig(figD, fullfile(fig_dir,'D_psd_ratio.png'), dpi);

% ===================== E) Band-limited RMSE bar with labels ===============
rmse_vec = [bands.rmse_low, bands.rmse_mid, bands.rmse_high];
share = (rmse_vec.^2) / max(eps, sum(rmse_vec.^2));   % energy shares (approx.)
figE = figure('Name','E_band_rmse','Color','w','Units','pixels'); figE.Position(3:4)=[720 420];
axE = axes(figE); hold(axE,'on'); set(axE,'LineWidth',axesLW);
bh = bar(axE, rmse_vec);
bh.FaceColor = 'flat';
bh.CData = [C.blue; C.purp; C.red];
set(axE,'XTickLabel',{'Low','Mid','High'});
ylabel(axE,'RMSE (m)'); title(axE,'Band-limited RMSE');
grid(axE,'on'); if useMinorGrid, grid(axE,'minor'); end
% Value + share labels
for i=1:numel(rmse_vec)
    text(i, rmse_vec(i), sprintf('%.4g m\n(%.0f%%)', rmse_vec(i), 100*share(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
end
save_fig(figE, fullfile(fig_dir,'E_band_rmse.png'), dpi);

% ===================== F) Coherence (optional) ============================
if have_mscohere && ~isempty(Cxy)
    kcoh = find(fcoh>0,1,'first'); fco = fcoh(kcoh:end); Cxy2 = Cxy(kcoh:end);
    figF = figure('Name','F_coherence','Color','w','Units','pixels'); figF.Position(3:4)=[980 420];
    tloF = tiledlayout(figF,1,1,'Padding','compact','TileSpacing','compact');
    axF = nexttile(tloF); hold(axF,'on'); set(axF,'LineWidth',axesLW,'XScale','log');
    plot(axF, fco, Cxy2, 'Color', C.green);
    yline(axF, 0.8, ':', 'Color', C.gray, 'DisplayName','0.8');
    yline(axF, 0.5, ':', 'Color', C.gray, 'DisplayName','0.5');
    grid(axF,'on'); if useMinorGrid, grid(axF,'minor'); end
    xlabel(axF,'Spatial frequency (cycles/m)'); ylabel(axF,'Coherence');
    ylim(axF, [0 1]); title(axF,'Magnitude-squared coherence (True vs Recon)');
    drawnow; yl = ylim(axF);
    shade_bands(axF, fco, yl, bands, nyq_local, 0.03);
    save_fig(figF, fullfile(fig_dir,'F_coherence.png'), dpi);
end

% ===================== G) Residual histogram ==============================
figG = figure('Name','G_residual_hist','Color','w','Units','pixels'); figG.Position(3:4)=[720 420];
axG = axes(figG); hold(axG,'on'); set(axG,'LineWidth',axesLW);
histogram(axG, align.residual, 50, 'Normalization','pdf','FaceColor',C.red, 'FaceAlpha',0.6, 'EdgeColor','none');
mu = mean(align.residual); sg = std(align.residual);
xline(axG, 0, ':', 'Color', C.gray, 'DisplayName','0');
xline(axG, mu, '-', 'Color', C.black, 'DisplayName','mean');
xline(axG, mu+sg, ':', 'Color', C.gray, 'HandleVisibility','off');
xline(axG, mu-sg, ':', 'Color', C.gray, 'HandleVisibility','off');
xline(axG, mu+2*sg, ':', 'Color', C.gray, 'HandleVisibility','off');
xline(axG, mu-2*sg, ':', 'Color', C.gray, 'HandleVisibility','off');
grid(axG,'on'); if useMinorGrid, grid(axG,'minor'); end
xlabel(axG,'Residual (m)'); ylabel(axG,'PDF');
title(axG, sprintf('Residual distribution — mu = %.4g m, sigma = %.4g m', mu, sg));
save_fig(figG, fullfile(fig_dir,'G_residual_hist.png'), dpi);

disp('[profile_registration_and_metrics] Figures saved.');

%% =========================
%  [10] HELPERS
%  =========================
function p = first_existing(cands)
    p = ''; for i=1:numel(cands), if exist(cands{i},'file'), p=cands{i}; return; end, end
end

function val = get_first(S, names, defaultVal)
    val = [];
    for i=1:numel(names)
        if isfield(S, names{i}) && ~isempty(S.(names{i})), val = S.(names{i}); return; end
    end
    if isempty(val), val = defaultVal; end
end

function y = choose_first_field(S, names)
    for i=1:numel(names)
        if isfield(S, names{i}) && ~isempty(S.(names{i}))
            y = S.(names{i})(:); return;
        end
    end
    error('None of the fields found: %s', strjoin(names,', '));
end

function [Pxx, f] = welch_psd_safe(x, fs_spatial, nfft_max)
    x = x(:); x = detrend(x,'constant');
    N = numel(x);
    if N < 8, Pxx=[]; f=[]; return; end
    nperseg  = min(N, max(128, 2^floor(log2(max(32, floor(N/4))))));
    noverlap = floor(nperseg/2);
    nfft     = 2^nextpow2(max(nperseg, min(nfft_max, N)));
    [Pxx, f] = pwelch(x, hamming(nperseg), noverlap, nfft, fs_spatial, 'onesided');
end

function s = ssim_1d(x, y)
    x=x(:); y=y(:); N=min(numel(x),numel(y)); x=x(1:N); y=y(1:N);
    K1=0.01; K2=0.03; L=max([range(x), range(y), eps]); C1=(K1*L)^2; C2=(K2*L)^2;
    win_len = max(11, 2*floor(0.25*N/10)+1); g = gausswin(win_len); g=g/sum(g);
    mu_x = conv(x,g,'same'); mu_y = conv(y,g,'same');
    sx2 = conv(x.^2,g,'same') - mu_x.^2;
    sy2 = conv(y.^2,g,'same') - mu_y.^2;
    sxy = conv(x.*y,g,'same') - mu_x.*mu_y;
    num = (2*mu_x.*mu_y + C1) .* (2*sxy + C2);
    den = (mu_x.^2 + mu_y.^2 + C1) .* (sx2 + sy2 + C2);
    s_map = num ./ max(den, eps);
    s = mean(s_map(~isnan(s_map) & ~isinf(s_map)));
end

function [x_true, zr_true] = extract_truth(T)
% Robustly read ground-truth road profile from simulation_results.mat.
% Supports your current 'results' struct and older 'sim' or top-level vars.

    x_true = []; zr_true = [];

    % --- NEW: your current structure ---
    if isfield(T,'results')
        R = T.results;

        % Prefer explicit spatial road grid + spatial profile
        if isfield(R,'spatial_road') && isfield(R,'zr') && isfield(R.zr,'x') ...
                && ~isempty(R.spatial_road) && ~isempty(R.zr.x)
            x_true  = R.spatial_road(:);
            zr_true = R.zr.x(:);
            return;
        end

        % If spatial_road missing but zr.x exists, build x from alternatives
        if isfield(R,'zr') && isfield(R.zr,'x') && ~isempty(R.zr.x)
            zr_true = R.zr.x(:);

            % 1) Use x_m if same length
            if isfield(R,'x_m') && ~isempty(R.x_m) && numel(R.x_m)==numel(zr_true)
                x_true = R.x_m(:);
                return;
            end

            % 2) Use known dx or fs_spatial
            dx = [];
            if isfield(R,'zr') && isfield(R.zr,'dx') && ~isempty(R.zr.dx)
                dx = R.zr.dx;
            elseif isfield(R,'meta') && isfield(R.meta,'fs_spatial') && ~isempty(R.meta.fs_spatial)
                dx = 1 / max(R.meta.fs_spatial, eps);
            end
            if ~isempty(dx)
                x0 = 0;
                x_true = x0 + (0:numel(zr_true)-1).' * dx;
                return;
            end

            % 3) Last resort: derive from time & speed (if present)
            if isfield(R,'time') && ~isempty(R.time) && isfield(R,'meta') && isfield(R.meta,'V_kmh') && ~isempty(R.meta.V_kmh)
                V = R.meta.V_kmh / 3.6;
                t = R.time(:);
                % Align lengths if off-by-one
                if numel(t) ~= numel(zr_true)
                    n = min(numel(t), numel(zr_true));
                    t = t(1:n); zr_true = zr_true(1:n);
                end
                x_true = t * V;
                return;
            end
        end
    end

    % --- Legacy/top-level formats (kept for compatibility) ---
    if isfield(T,'x_true') && isfield(T,'zr_true')
        x_true  = T.x_true(:);
        zr_true = T.zr_true(:);
        return;
    end
    if isfield(T,'sim')
        sim = T.sim;
        if isfield(sim,'x_true') && isfield(sim,'zr_true')
            x_true  = sim.x_true(:);
            zr_true = sim.zr_true(:);
            return;
        end
        if isfield(sim,'t') && isfield(sim,'zr_true') && isfield(sim,'v')
            t = sim.t(:);
            v = sim.v(:); if isscalar(v), v = v*ones(size(t)); end
            x_true  = cumtrapz(t, v);
            zr_true = sim.zr_true(:);
            return;
        end
    end

    error('Could not resolve true profile fields in simulation_results.mat');
end

function s = truncate_path(p, maxlen)
    if numel(p)<=maxlen, s=p; else, s=['...' p(end-maxlen+4:end)]; end
end

function save_fig(fig, filepath, dpi)
% Clean export without transparency warnings
    try
        exportgraphics(fig, filepath, 'Resolution', dpi);
    catch
        [p,f,e] = fileparts(filepath);
        if isempty(e), e = '.png'; end
        switch lower(e)
            case '.png',  print(fig, fullfile(p,[f e]), '-dpng',  sprintf('-r%d',dpi));
            case '.tif',  print(fig, fullfile(p,[f e]), '-dtiff', sprintf('-r%d',dpi));
            case '.pdf'
                set(fig,'PaperPositionMode','auto'); 
                print(fig, fullfile(p,[f e]), '-dpdf','-vector');
            otherwise
                print(fig, fullfile(p,[f '.png']), '-dpng', sprintf('-r%d',dpi));
        end
    end
    try savefig(fig, [filepath(1:end-4) '.fig']); catch, end
end

function shade_bands(ax, f, ylims, bands, nyq, alpha)
% Shade low/mid/high spatial-frequency bands on an axes with x = f (cycles/m)
% alpha is optional (default 0.08)
    if nargin<6 || isempty(alpha), alpha = 0.08; end
    hold(ax,'on');
    Cb = [0.80 0.90 1.00];  % low (blueish)
    Cm = [0.92 0.85 0.96];  % mid (purple-ish)
    Ch = [1.00 0.88 0.88];  % high (reddish)

    xmin = f(1); xmax = f(end);
    bands_l = [0, min(bands.c_low.max, 0.999*nyq)];
    bands_m = [max(bands.c_mid.min,0), min(bands.c_mid.max, 0.999*nyq)];
    bands_h = [max(bands.c_high.min,0), 0.999*nyq];

    draw_patch(ax, max(xmin,bands_l(1)), min(xmax,bands_l(2)), ylims, Cb, alpha, 'Low');
    draw_patch(ax, max(xmin,bands_m(1)), min(xmax,bands_m(2)), ylims, Cm, alpha, 'Mid');
    draw_patch(ax, max(xmin,bands_h(1)), min(xmax,bands_h(2)), ylims, Ch, alpha, 'High');

    % Band separators with labels
    if isfield(bands,'c_low') && isfield(bands.c_low,'max')
        xline(ax, bands.c_low.max, ':', 'Low/Mid', 'Color',[0.4 0.4 0.4], 'LabelVerticalAlignment','bottom');
    end
    if isfield(bands,'c_mid') && isfield(bands.c_mid,'max')
        xline(ax, bands.c_mid.max, ':', 'Mid/High', 'Color',[0.4 0.4 0.4], 'LabelVerticalAlignment','bottom');
    end
end

function draw_patch(ax, x1, x2, ylims, color, alpha, ~)
    if ~isfinite(x1) || ~isfinite(x2) || x2<=x1, return; end
    fx = [x1 x2 x2 x1]; fy = [ylims(1) ylims(1) ylims(2) ylims(2)];
    patch('Parent',ax,'XData',fx,'YData',fy,'FaceColor',color,'FaceAlpha',alpha, ...
          'EdgeColor','none','HandleVisibility','off');
end

function yL = symmetric_limits(y)
% Symmetric limits around zero using robust max (98th percentile)
    y = y(isfinite(y));
    if isempty(y), yL = [-1 1]; return; end
    m = prctile(abs(y),98); 
    m = max(m, max(abs(y))*0.5); % avoid too-tight if heavy tails
    yL = 1.05*[-m m];
end