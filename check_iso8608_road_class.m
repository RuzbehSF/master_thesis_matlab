
%% Validating Road Class using ISO 8608
% -------------------------------------------------------------------------
% Verifies the ISO 8608 road class from the saved results of your
% quarter-car simulation by estimating spatial PSD Gd(n), fitting waviness
% w and Gd(n0) at n0=0.1 cycles/m, and comparing to class thresholds.
% Produces a log-log PSD plot with fit band, confidence band, class swatch,
% residuals subplot, and saves a PNG next to the MAT file.
% -------------------------------------------------------------------------

clear; clc; close all;

%% 0) Where to load results
results_file = fullfile('00_Outputs','01_Simulation','SimulationData','simulation_results.mat');
assert(exist(results_file,'file')==2, ...
    'Could not find results file at: %s. Run your main script first.', results_file);

S = load(results_file);
results = S.results;

% Pull what we need
x   = results.spatial_road(:);         % m
if isfield(results.zr,'base_x')
    zr = results.zr.base_x(:);
else
    zr = results.zr.x(:);
end
dx0 = results.zr.dx;                   % m (nominal spacing)
fsx_nom = 1/dx0;                       % samples per meter (cycles/m sampling rate)
declared_class = upper(string(results.meta.road_class));

% Quick uniform spacing sanity check
dx = diff(x);
dx_err = max(abs(dx - median(dx)));
if dx_err > 1e-6
    warning('Spatial step is not perfectly uniform (max dev = %.3g m). Using median dx.', dx_err);
end
dx  = median(dx);
fsx = results.meta.fs_spatial;         % samples per meter (from simulation metadata)

% Detrend to remove any weak drift (keeps ISO band shape)
zr = detrend(zr,'linear');

%% 1) Welch PSD in the spatial domain (units: m^3 per (cycles/m))
% Choose a window length ~ 8–16k samples if available; otherwise adapt.
N = numel(zr);
nfft   = 2^nextpow2(min(N, 2^16));
win_len= min(max(2^nextpow2(round(N/8)), 1024), N);
win    = hann(win_len,'periodic');
nover  = round(0.5*win_len);

% Spatial PSD: pwelch with fs argument as "samples per meter" (cycles/m)
[Gd_hat, n_cycperm] = pwelch(zr, win, nover, nfft, fsx, 'onesided');

% Keep a reasonable ISO band for fitting (typical: 0.02–3 cycles/m)
fit_band = [0.011, 2.828];  % ISO 8608 reporting band used in your generator
m_fit  = n_cycperm >= fit_band(1) & n_cycperm <= fit_band(2);
n_fit  = n_cycperm(m_fit);
Gd_fit = Gd_hat(m_fit);

% Guard against empties
assert(~isempty(n_fit) && numel(n_fit)>10, ...
    'Not enough PSD points in fit band [%.3g, %.3g] cycles/m.', fit_band);

%% 2) Fit ISO line: log10(Gd) = alpha + beta * log10(n)
n0   = results.iso8608.n0_cycpm;  % cycles/m from simulation metadata
xlog = log10(n_fit);
ylog = log10(Gd_fit);
P    = polyfit(xlog, ylog, 1);    % [beta, alpha]
beta  = P(1);
alpha = P(2);

w_est   = -beta;
Gd0_est = 10^(alpha + beta*log10(n0));

%% 3) Class decision using ISO-8608-style centers and geometric boundaries
class_names   = ["A","B","C","D","E","F","G","H"];
class_centers = results.iso8608.center_map(:).'; % m^3 (A..H)
lower_bounds  = results.iso8608.bounds_lo(:).';  % geometric lower bounds
upper_bounds  = results.iso8608.bounds_hi(:).';  % geometric upper bounds

% Locate class for Gd0_est
class_idx = find(Gd0_est >= lower_bounds & Gd0_est < upper_bounds, 1, 'first');
if isempty(class_idx), class_idx = numel(class_names); end
estimated_class = class_names(class_idx);

%% 4) Pass/Fail gates and textual report
w_tol = 0.25;  % acceptable deviation around 2.0
class_match = (estimated_class == declared_class);
w_ok = abs(w_est - 2.0) <= w_tol;

fprintf('\n=== ISO 8608 Check ===\n');
fprintf('Declared class:    %s\n', declared_class);
fprintf('Estimated class:   %s\n', estimated_class);
fprintf('Estimated w:       %.3f  (target ~ 2.0; tol ±%.2f) %s\n', ...
    w_est, w_tol, ternary(w_ok,'OK','(check)'));
fprintf('Estimated Gd(n0):  %.3e  m^3  at n0 = %.2f cycles/m\n', Gd0_est, n0);
fprintf('Class center Gd0:  %.3e  m^3  for class %s\n', class_centers(class_idx), estimated_class);

if class_match && w_ok
    fprintf('RESULT: PASS ✅ — PSD level and waviness consistent with class %s.\n', declared_class);
elseif class_match && ~w_ok
    fprintf('RESULT: WARN ⚠️ — Gd(n0) matches class %s, but waviness w=%.2f deviates from 2.0±%.2f.\n', ...
        declared_class, w_est, w_tol);
elseif ~class_match && w_ok
    fprintf('RESULT: MISMATCH ❌ — Waviness OK, but class likely %s (declared %s).\n', ...
        estimated_class, declared_class);
else
    fprintf('RESULT: MISMATCH ❌ — Both class and waviness deviate; revisit road synthesis.\n');
end

%% 5) Plot  (improved visuals + legend fix)
figure('Name','ISO 8608 PSD Check','Position',[80 80 1080 720]); 
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ---------- (a) PSD on log-log with fit, fit-band, CI, class swatch ----------
ax1 = nexttile(1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');

% 95% CI for Welch PSD (approx) ---------------------------------------------
K_est = max(1, floor((N - nover)/(max(1,win_len - nover))));  % num segments (rough)
dof = 2*K_est;                                                % chi-square DOF
show_CI = (dof >= 4) && exist('chi2inv','file')==2;
alpha_CI = 0.05;

hCI = gobjects(1,1);  % predeclare handle
if show_CI
    L = dof .* Gd_hat ./ chi2inv(1 - alpha_CI/2, dof);
    U = dof .* Gd_hat ./ chi2inv(alpha_CI/2, dof);
    fill_x = [n_cycperm; flipud(n_cycperm)];
    fill_y = [L; flipud(U)];
    hCI = fill(fill_x, fill_y, [0.8 0.85 1.0], ...
        'EdgeColor','none', 'FaceAlpha',0.35, 'DisplayName','95% CI (Welch)');
end

% Fit band shading -----------------------------------------------------------
ymin_plot = min(Gd_hat(Gd_hat>0))*0.7;
ymax_plot = max(Gd_hat)*1.4;
patch([fit_band(1) fit_band(2) fit_band(2) fit_band(1)], ...
      [ymin_plot ymax_plot ymax_plot ymin_plot], ...
      [0.93 0.93 0.93], 'EdgeColor','none','FaceAlpha',0.5, 'DisplayName','Fit band');

% PSD curve ------------------------------------------------------------------
hPSD = loglog(n_cycperm, Gd_hat, 'LineWidth', 1.4, 'DisplayName','G_d(n)');

% Fitted line over fit band --------------------------------------------------
n_plot = logspace(log10(fit_band(1)), log10(fit_band(2)), 300);
Gd_fitline = 10.^(alpha + beta*log10(n_plot));
hFit = loglog(n_plot, Gd_fitline, '--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Fit: w = %.2f', w_est));

% Vertical line at n0 --------------------------------------------------------
xline(n0, ':k', 'n_0 = 0.1 cyc/m', 'LabelVerticalAlignment','bottom', ...
    'LabelOrientation','horizontal');

% Marker at n0 with estimated level -----------------------------------------
Gd_at_n0 = 10^(alpha + beta*log10(n0));
plot(n0, Gd_at_n0, 'o', 'MarkerSize', 7, 'LineWidth',1.2, ...
     'DisplayName','G_d(n_0) (est)');

% Reference slope guide (w = 2) through (n0, Gd_at_n0) ----------------------
guide_w = 2.0;
guide_alpha = log10(Gd_at_n0) + guide_w*log10(n0); % passes through (n0,Gd_at_n0)
Gd_guide = 10.^(guide_alpha + (-guide_w)*log10(n_plot));
loglog(n_plot, Gd_guide, ':', 'LineWidth', 1.4, 'DisplayName','Ref slope w = 2');

% ISO class swatch around n0 (visual comparator) -----------------------------
strip_w = [n0/1.6, n0*1.6];  % narrow vertical strip
for k = 1:numel(class_centers)
    g0 = class_centers(k);
    plot(strip_w, [g0 g0], '-', 'LineWidth', 1.0, 'Color',[0.45 0.45 0.45]);
end
% Highlight the estimated class band
lb = lower_bounds(max(1, class_idx));  
ub = upper_bounds(min(numel(upper_bounds), class_idx));
patch([strip_w(1) strip_w(2) strip_w(2) strip_w(1)], [lb lb ub ub], ...
      [0.95 0.90 0.80], 'EdgeColor',[0.4 0.3 0.1], 'LineWidth',0.8, 'FaceAlpha',0.35, ...
      'DisplayName',sprintf('Class %s band @ n_0', estimated_class));

% Cosmetics ------------------------------------------------------------------
xlabel('Spatial frequency n (cycles/m)');
ylabel('G_d(n)  [m^3]');
title(sprintf('ISO 8608 PSD — Declared %s | Estimated %s | w = %.2f', ...
    declared_class, estimated_class, w_est));

xlim([max(min(n_cycperm(n_cycperm>0))*0.9, 1e-3), min(max(n_cycperm)*1.1, 10)]);
ylim([ymin_plot, max(Gd_hat)*1.5]);

% Legend: include CI handle only if plotted (fixes identical-block warning)
leg_handles = [hPSD hFit];
if show_CI && isgraphics(hCI)
    leg_handles = [hCI leg_handles];
end
legend(leg_handles, 'Location','southwest');

set(ax1,'FontName','Helvetica','FontSize',11,'LineWidth',0.9);

% ---------- (b) Residuals subplot (how power-law fits across band) ----------
ax2 = nexttile(2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');

% Residuals in dB (log10 space): 10*log10(G_hat / G_fit)
G_fit_band = 10.^(alpha + beta*log10(n_fit));
res_dB = 10*log10(Gd_fit ./ G_fit_band);

plot(n_fit, res_dB, '-', 'LineWidth', 1.3);
yline(0, ':k','Baseline');
yline(+3, ':r'); 
yline(-3, ':r');  % ±3 dB visual tolerance

xlabel('Spatial frequency n (cycles/m)');
ylabel('Residual (dB)');
title('Fit residuals over band');
set(ax2,'XScale','log');

xlim([fit_band(1) fit_band(2)]);
ylim(max(abs(res_dB))*[-1.2 1.2] + [-0.5 0.5]);  % small padding
set(ax2,'FontName','Helvetica','FontSize',11,'LineWidth',0.9);

% ---------- Optional: save a high-res copy next to results ------------------
try
    out_png = fullfile('00_Outputs','01_Simulation','SimulationData','iso8608_psd_check.png');
    exportgraphics(gcf, out_png, 'Resolution', 200);
    fprintf('Saved PSD check figure: %s\n', out_png);
catch
    % no-op
end

%% 6) Table summary
tbl = table( ...
    string(declared_class), string(estimated_class), w_est, Gd0_est, ...
    'VariableNames', {'DeclaredClass','EstimatedClass','Waviness_w','Gd_n0_m3'});
disp(' ');
disp(tbl);

%% Helper: inline ternary
function out = ternary(cond,a,b)
if cond
    out = a;
else
    out = b;
end
end