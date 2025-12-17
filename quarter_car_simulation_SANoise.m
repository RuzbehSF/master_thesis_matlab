
%% 1) Simulate Quarter-Car Vehicle Dynamics Over a Generated Road
% -------------------------------------------------------------------------
% PURPOSE:
% This script is the starting point of the analysis pipeline. It simulates
% the vertical dynamic response of a 2-DOF (two-degree-of-freedom)
% quarter-car model as it travels over a procedurally generated road
% profile. Its primary goal is to create a realistic dataset, including
% vehicle acceleration and the 'ground truth' road profile, which
% subsequent scripts will use.
%
% CORE LOGIC (How it works):
% 1.  Defines all physical constants for the vehicle (using 'Golden Car'
%     benchmark parameters) and settings for the simulation (e.g., speed,
%     road length, sampling rate).
%
% 2.  Procedurally generates a road surface by combining a stochastic
%     background roughness (based on ISO 8608 classes) with a library of
%     specific damage types (e.g., potholes, bumps) selected via a
%     'damage_scenario' setting.
%
% 3.  Represents the quarter-car's physics as a linear state-space model
%     (ẋ = Ax + Bu), where the state vector 'x' includes the positions
%     and velocities of the sprung and unsprung masses.
%
% 4.  Solves the state-space equations over time using a numerical solver
%     (e.g., lsim). This computes the vehicle's motion (displacements,
%     velocities, accelerations) at each time step in response to the
%     road input.
%
% 5.  Organizes all inputs, parameters, and time-history results into a
%     single structured .mat file, which serves as the ground-truth
%     dataset for the rest of the pipeline.
%
% INPUTS:
% - None. All parameters are defined within the script.
%
% OUTPUTS:
% - '01_Quarter Car Modeling/simulation_results.mat': A MAT-file
%   containing all simulation parameters, time-domain results, and road
%   profile data.
% - Figures visualizing the road profile and vehicle response.
% -------------------------------------------------------------------------

close all; clear; clc;

%% 0) Hyperparameters and Simulation Settings
%%% Vehicle Preset (per-corner parameters via getVehiclePreset)

% Available presets (quick guidance):
% 'Golden_Car'   : IRI benchmark. No tire damping (ct = 0). Good for apples-to-apples comparisons.
% 'Sport_Sedan'  : Lower sprung mass, stiffer ks, higher cs, stiffer tire → taut, higher natural freqs.
% 'SUV'          : Higher sprung mass; comfort-leaning ks/cs; slightly softer tire than sport sedan.
% 'Light_Truck'  : Higher masses and ks to support payload; moderate-high damping; heavier unsprung mass.
% 'EV_Compact'   : Heavier sprung mass (battery); balanced ks/cs; slightly stiffer tire; small tire damping.
% 'Luxury_Sedan' : Comfort tuned—softer ks and cs; tire stiffness similar to baseline; plush primary ride.
% 'Offroad_4x4'  : Larger unsprung mass; softer ks for travel, higher cs for control; compliant tire.
vehicle_preset  = 'Golden_Car';

% Optional ad-hoc overrides (leave empty to use preset as-is), e.g.:
% vehicle_override = struct('Suspension_Damping', 8.0e3, 'Tire_Damping', 200);
vehicle_override = struct();

% Load vehicle parameters from preset
[VH, ~] = getVehiclePreset(vehicle_preset);

% Assign required variables (names preserved for the rest of the script)
Sprung_Mass           = VH.Sprung_Mass;          % kg  (per corner sprung mass)
Unsprung_Mass         = VH.Unsprung_Mass;        % kg  (per wheel/tire/axle)
Suspension_Stiffness  = VH.Suspension_Stiffness; % N/m (spring)
Suspension_Damping    = VH.Suspension_Damping;   % Ns/m (damper)
Tire_Stiffness        = VH.Tire_Stiffness;       % N/m (tire)
Tire_Damping          = VH.Tire_Damping;         % Ns/m

% Apply overrides if provided (only fields present will replace the above)
if ~isempty(vehicle_override)
    if isfield(vehicle_override, 'Sprung_Mass')
        Sprung_Mass = vehicle_override.Sprung_Mass;
    end
    if isfield(vehicle_override, 'Unsprung_Mass')
        Unsprung_Mass = vehicle_override.Unsprung_Mass;
    end
    if isfield(vehicle_override, 'Suspension_Stiffness')
        Suspension_Stiffness = vehicle_override.Suspension_Stiffness;
    end
    if isfield(vehicle_override, 'Suspension_Damping')
        Suspension_Damping = vehicle_override.Suspension_Damping;
    end
    if isfield(vehicle_override, 'Tire_Stiffness')
        Tire_Stiffness = vehicle_override.Tire_Stiffness;
    end
    if isfield(vehicle_override, 'Tire_Damping')
        Tire_Damping = vehicle_override.Tire_Damping;
    end
end

% total vehicle mass from quarter-car values:
Total_Mass = 4*(Sprung_Mass + Unsprung_Mass);   % ≈ vehicle curb mass (kg)

%%% Road & Simulation
Road_Length   = 2000;       % m
Vehicle_Speed = 80;         % km/h
Sampling_Freq = 200;        % Hz

%%% Smartphone Accelerometer Noise Modeling

enable_noise = false;

% 1. Broadband White Noise (Noise Density)
% (Typical smartphone MEMS: 150-400 ug/sqrt(Hz))
% g/sqrt(Hz) (Bosch BMI160, [Lăpădat et al., 2021])
noise_density_g_sqrt_hz = 180e-6;

% 2. Bias Drift (Random Walk)
% (This models the slow-wandering offset)
% (standard deviation of random walk step, calibrated to Fig A1 of [Lăpădat et al., 2021])a
bias_drift_std_g_per_step = 1.8e-6;       % g

% 3. Static Bias Offset
%    (A constant offset from zero, e.g., 20 mg)
bias_offset_g = 0.02;                     % g

%%% Road Roughness Class (ISO 8608)
% Classes are defined by the road elevation PSD level Gd(n0) at n0 = 0.1 cycles/m (units: m^3).
% Reference mapping (typical values):
%   A : 16e-6     → Very smooth / new pavement
%   B : 64e-6     → Very good
%   C : 256e-6    → Good
%   D : 1024e-6   → Average
%   E : 4096e-6   → Poor
%   F : 16384e-6  → Very poor
%   G : 65536e-6  → Extremely poor
%   H : 262144e-6 → Severe / off-road
road_class = 'C';

%%% Damage scenario selection:
% Scenarios (exact names):
% 'Smooth'                       (Base ISO roughness only, no damage features)
% 'Bumps_Only'                   (4 bumps of increasing severity: low → medium → high → extreme)
% 'Potholes_Only'                (4 potholes of increasing severity: low → medium → high → extreme)
% 'Rutting_Only'                 (4 rut sections of increasing severity: low → medium → high → extreme)
% 'Crack_Joints_Only'            (4 individual crack/joint features of increasing severity)
% 'Block_Cracking_Only'          (4 road sections with periodic block cracking patterns)
% 'Fatigue_Crocodile_Only'       (4 road sections with fatigue/alligator cracking fields)
% 'Reflection_Cracking_Only'     (4 road sections with periodic reflection cracking)
% 'Corrugation_Shoving_Only'     (4 road sections with corrugation ripples + shoving steps)
% 'Depression_Only'              (4 local depression features of increasing severity)
% 'Swelling_Only'                (4 local swelling features of increasing severity)
% 'Edge_Fatigue_Only'            (4 road sections with edge fatigue damage patterns)

%%%% New Scenarios:
% 'Highway_Smooth'               (Realistic highway conditions with expansion joints and utility cuts)
% 'City_Urban'                   (Urban driving with speed bumps, potholes, and utility patches)
% 'Progressive_Roughness'        (Road class A → H progression, no damage features)

%%%% Severity-Based Scenarios:
% 'Low_Severity'                 (All damage types at low severity level)
% 'Medium_Severity'              (All damage types at medium severity level)
% 'High_Severity'                (All damage types at high severity level)
% 'Extreme_Severity'             (All damage types at extreme severity level)
% 'All_Damages_Mixed_Severity'   (All damage types with mixed severity levels)

damage_scenario_selection = 'Progressive_Roughness';

%%% ODE Solver Options
ode_RelTol = 1e-3;
ode_AbsTol = 1e-4;

% Solver method:
% 'ode45'
% 'ode15s'
% 'ss-cont'
% 'ss-disc'
% 'ss-disc-manual'
solver_method = 'ss-disc';

save_results = true;

%% 1) Derived sim params
V  = Vehicle_Speed / 3.6;                                % m/s
dt = 1 / Sampling_Freq;                                  % s
time_road_profile = 0 : dt : (Road_Length / V);

%% 2) Quick modal sanity check (quarter-car linearization)
ms = Sprung_Mass;
mus = Unsprung_Mass;
ks = Suspension_Stiffness;
cs = Suspension_Damping;
kt = Tire_Stiffness;
ct = Tire_Damping;
A_chk = [0 1 0 0; -ks/ms -cs/ms ks/ms cs/ms; 0 0 0 1; ks/mus cs/mus -(ks+kt)/mus -(cs+ct)/mus];
ev = eig(A_chk);
wn = abs(imag(ev));
zeta = -real(ev)./sqrt(real(ev).^2+imag(ev).^2+eps); f_hz = wn/(2*pi);
[fsorted,idxm] = sort(f_hz);
% Collapse conjugate duplicates to unique modes (tolerance 0.01 Hz)
[~, ia] = unique(round(fsorted,2),'stable');
fu = fsorted(ia);
zu = zeta(idxm(ia));
if numel(fu) >= 2
    fprintf('Quarter-car modes (Hz,zeta): %.2f(%.2f), %.2f(%.2f)\n', fu(1), zu(1), fu(2), zu(2));
    fmax = fu(2);
elseif isscalar(fu)
    fprintf('Quarter-car mode (Hz,zeta): %.2f(%.2f)\n', fu(1), zu(1));
    fmax = fu(1);
else
    fmax = 0;
end
if Sampling_Freq < 10*fmax && fmax > 0
    warning('Sampling_Freq = %.1f Hz may be low vs highest mode = %.2f Hz. Consider fs >= %.1f Hz.', ...
        Sampling_Freq, fmax, 10*fmax);
end

%% 3) Road profile (roughness + damages)
[spatial_road, z_r_x, time_sim, z_r_t, dz_r_t, z_r_x_base] = ...
    generateRoadProfile(Road_Length, V, dt, time_road_profile, road_class, damage_scenario_selection);

%% 4) Build road inputs for quarter-car model
[U, wheel_labels] = buildRoadInputs(spatial_road, z_r_x, time_sim, V);

%% 5) Build quarter-car model and simulate
params = packParams(Sprung_Mass, Unsprung_Mass, Suspension_Stiffness, Suspension_Damping, ...
                    Tire_Stiffness, Tire_Damping);

[A,B] = quarter_AB(params);
X0 = [0;0;0;0];
C = eye(4);
D = zeros(4,size(B,2));

switch lower(solver_method)
    case 'ode45'
        f = @(t,y) quarterODE(t,y,U,time_sim,params);
        [t_sol, Y_sol] = ode45(@(t,y) f(t,y), time_sim, X0, odeset('RelTol',ode_RelTol,'AbsTol',ode_AbsTol));
    case 'ode15s'
        f = @(t,y) quarterODE(t,y,U,time_sim,params);
        [t_sol, Y_sol] = ode15s(@(t,y) f(t,y), time_sim, X0, odeset('RelTol',ode_RelTol,'AbsTol',ode_AbsTol));
    case 'ss-cont'
        sys = ss(A,B,C,D); [Y_sol, t_sol] = lsim(sys, U, time_sim, X0);
    case 'ss-disc'
        sysc = ss(A,B,C,D); sysd = c2d(sysc, 1/Sampling_Freq, 'zoh');
        [Y_sol, t_sol] = lsim(sysd, U, time_sim, X0);
    case 'ss-disc-manual'
        sysc = ss(A,B,C,D); sysd = c2d(sysc, 1/Sampling_Freq, 'zoh');
        Ad = sysd.A;
        Bd = sysd.B;
        N = numel(time_sim);
        X = zeros(size(A,1), N);
        X(:,1) = X0(:);
        for k=1:N-1
            X(:,k+1) = Ad*X(:,k) + Bd*U(k,:)';
        end
        Y_sol = X.';
        t_sol = time_sim(:);
    otherwise
        error('Unknown solver_method: %s', solver_method);
end

%% 6) Body outputs (heave & accel) for quarter-car
[Body_Heave_m, Body_Accel_clean] = bodyOutputs(params, Y_sol, U, t_sol);

% Keep legacy variable names for your plots/save:
Sprung_Mass_Disp      = Body_Heave_m;       % m (body CG heave or ys in quarter)
sprung_mass_accel_raw = Body_Accel_clean;   % m/s^2

% Map states to animation-friendly variables
Unsprung_Mass_Disp = Y_sol(:,3);            % yus (m)

%% 7) Add realistic accelerometer noise (toggle-aware)
N = numel(sprung_mass_accel_raw);
g0 = 9.806;

if enable_noise
    % 1) Broadband white noise from noise density (g/√Hz → m/s²)
    nyquist_freq = Sampling_Freq / 2;
    noise_std_g   = noise_density_g_sqrt_hz * sqrt(nyquist_freq);
    noise_std_mps2 = noise_std_g * g0;
    white_noise   = noise_std_mps2 * randn(N, 1);

    % 2) Bias drift (random walk)
    bias_drift_std_mps2 = bias_drift_std_g_per_step * g0;
    bias_drift          = cumsum(bias_drift_std_mps2 * randn(N, 1));

    % 3) Static bias offset
    bias_offset_mps2 = bias_offset_g * g0;

    % 4) Sum all components
    sprung_mass_accel = sprung_mass_accel_raw + white_noise + bias_drift + bias_offset_mps2;
else
    % Noise disabled → keep the measured channel identical to the clean one
    white_noise       = zeros(N,1);
    bias_drift        = zeros(N,1);
    bias_offset_mps2  = 0;
    sprung_mass_accel = sprung_mass_accel_raw;
end

%% 8) Visualization
g0 = 9.806;
figName = sprintf('Vehicle: QUARTER | Scenario: %s', strrep(damage_scenario_selection,'_',' '));

% Create folder for saving plots
save_folder = fullfile('00_Outputs', '01_Simulation', 'SimulationData');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

figure('Name', figName, 'Position', [100, 70, 980, 700]);
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

% Road (distance)
ax1 = nexttile(1);
plot(spatial_road, z_r_x*1000,'LineWidth',1);
if strcmpi(damage_scenario_selection,'progressive_roughness')
    class_str = 'Classes A→H';
else
    class_str = sprintf('Class %s', road_class);
end
title(sprintf('Road Profile (%s): %s', class_str, strrep(damage_scenario_selection,'_',' ')));
xlabel('Distance (m)'); ylabel('Elevation (mm)');
grid on; grid minor; box on;
legend('Road','Location','best');
% Dynamic y-axis range based on data
road_range = max(abs(z_r_x*1000));
span = max(road_range*1.1, 1e-6);
ylim([-span, span]);

% Body heave (time)
ax2 = nexttile(2);
plot(t_sol, Sprung_Mass_Disp*1000,'LineWidth',1);
title('Body Heave (CG)');
xlabel('Time (s)');
ylabel('Displacement (mm)');
grid on; grid minor; box on;
% Dynamic y-axis range based on data
heave_range = max(abs(Sprung_Mass_Disp*1000));
span = max(heave_range*1.1, 1e-6);
ylim([-span, span]);

% Noisy accel
ax3 = nexttile(3);
yyaxis left; 
h1=plot(t_sol, sprung_mass_accel,'-','LineWidth',1);
ylabel('m/s^2');
hold on;
yline(0,'k-','LineWidth',0.5);
yyaxis right;
h2 = plot(t_sol, sprung_mass_accel/g0,'--','LineWidth',0.9);
ylabel('g');
title('Body Vertical Acceleration (Noisy)');
xlabel('Time (s)');
grid on; grid minor; box on;
legend([h1 h2],{'Accel (m/s^2)','Accel (g)'},'Location','best');
% Dynamic y-axis range based on data
accel_range = max(abs(sprung_mass_accel));
yyaxis left;
spanL = max(accel_range * 1.1 , 1e-6);
ylim([-spanL, spanL]);
yyaxis right;
spanR = max(accel_range * 1.1/g0 , 1e-9);
ylim([-spanR, spanR]);

% Clean accel
ax4 = nexttile(4);
plot(t_sol, sprung_mass_accel_raw,'k-','LineWidth',1);
hold on;
yline(0,'k-','LineWidth',0.5);
title('Body Vertical Acceleration (Clean)');
xlabel('Time (s)');
ylabel('m/s^2');
grid on; grid minor; box on;
% Dynamic y-axis range based on data
clean_accel_range = max(abs(sprung_mass_accel_raw));
span = max(clean_accel_range*1.1, 1e-6);
ylim([-span, span]);

linkaxes([ax2,ax3,ax4],'x');

% Save main figure
if save_results
    main_fig_path = fullfile(save_folder, 'simulation_results.png');
    saveas(gcf, main_fig_path);
    disp(['Saved: ' main_fig_path]);
end

%% 8.1) Additional Analysis Figures
% % Frequency Domain Analysis
% figure('Name', [figName ' - Frequency Analysis'], 'Position', [150, 100, 1000, 600]);
% tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% 
% % Power Spectral Density
% ax5 = nexttile(1);
% [psd_accel, f] = pwelch(sprung_mass_accel_raw, [], [], [], Sampling_Freq);
% semilogy(f, psd_accel, 'b-', 'LineWidth', 1.5);
% title('Power Spectral Density (Clean Acceleration)','FontSize',12,'FontWeight','bold');
% xlabel('Frequency (Hz)','FontSize',11); ylabel('PSD (m²/s⁴/Hz)','FontSize',11);
% grid on; grid minor; box on;
% xlim([0, 50]); % Focus on 0-50 Hz range
% % Add quarter-car natural frequencies
% wn1 = sqrt((Suspension_Stiffness/Sprung_Mass) + (Tire_Stiffness/Unsprung_Mass));
% wn2 = sqrt(Suspension_Stiffness/Sprung_Mass);
% xline(wn1/(2*pi), 'r--', sprintf('Body Mode: %.1f Hz', wn1/(2*pi)), 'LineWidth', 1);
% xline(wn2/(2*pi), 'g--', sprintf('Wheel Mode: %.1f Hz', wn2/(2*pi)), 'LineWidth', 1);
% 
% % Road Profile Statistics
% ax6 = nexttile(2);
% histogram(z_r_x*1000, 50, 'FaceColor', 'cyan', 'EdgeColor', 'black', 'FaceAlpha', 0.7);
% title('Road Profile Elevation Distribution','FontSize',12,'FontWeight','bold');
% xlabel('Elevation (mm)','FontSize',11); ylabel('Frequency','FontSize',11);
% grid on; grid minor; box on;
% % Add statistics
% mean_elev = mean(z_r_x*1000);
% std_elev = std(z_r_x*1000);
% rms_elev = rms(z_r_x*1000);
% text(0.02, 0.98, sprintf('Mean: %.2f mm\nStd: %.2f mm\nRMS: %.2f mm', mean_elev, std_elev, rms_elev), ...
%      'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','white');
% 
% % Acceleration vs Displacement (Phase Plot)
% ax7 = nexttile(3);
% plot(Sprung_Mass_Disp*1000, sprung_mass_accel_raw, 'b-', 'LineWidth', 1);
% title('Phase Plot: Acceleration vs Displacement','FontSize',12,'FontWeight','bold');
% xlabel('Displacement (mm)','FontSize',11); ylabel('Acceleration (m/s²)','FontSize',11);
% grid on; grid minor; box on;
% 
% % Cumulative RMS
% ax8 = nexttile(4);
% cumulative_rms = sqrt(cumsum(sprung_mass_accel_raw.^2) ./ (1:length(sprung_mass_accel_raw))');
% plot(t_sol, cumulative_rms, 'r-', 'LineWidth', 1.5);
% title('Cumulative RMS Acceleration','FontSize',12,'FontWeight','bold');
% xlabel('Time (s)','FontSize',11); ylabel('RMS Acceleration (m/s²)','FontSize',11);
% grid on; grid minor; box on;
% 
% % Save frequency analysis figure
% if save_results
%     freq_fig_path = fullfile(save_folder, 'simulation_frequency_analysis.png');
%     saveas(gcf, freq_fig_path);
%     disp(['Saved: ' freq_fig_path]);
% end

%% 8.2) Performance Metrics Figure
% figure('Name', [figName ' - Performance Metrics'], 'Position', [200, 150, 1000, 600]);
% tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% 
% % Road Profile Gradient
% ax9 = nexttile(1);
% road_gradient = gradient(z_r_x*1000, spatial_road);
% plot(spatial_road, road_gradient, 'm-', 'LineWidth', 1.5);
% title('Road Profile Gradient','FontSize',12,'FontWeight','bold');
% xlabel('Distance (m)','FontSize',11); ylabel('Gradient (mm/m)','FontSize',11);
% grid on; grid minor; box on;
% % Dynamic y-axis range based on data
% gradient_range = max(abs(road_gradient));
% ylim([-gradient_range*1.1, gradient_range*1.1]);
% 
% % Acceleration Histogram
% ax10 = nexttile(2);
% histogram(sprung_mass_accel_raw, 50, 'FaceColor', 'yellow', 'EdgeColor', 'black', 'FaceAlpha', 0.7);
% title('Acceleration Distribution','FontSize',12,'FontWeight','bold');
% xlabel('Acceleration (m/s²)','FontSize',11); ylabel('Frequency','FontSize',11);
% grid on; grid minor; box on;
% % Add normal distribution overlay
% hold on;
% x_norm = linspace(min(sprung_mass_accel_raw), max(sprung_mass_accel_raw), 100);
% y_norm = normpdf(x_norm, mean(sprung_mass_accel_raw), std(sprung_mass_accel_raw));
% y_norm = y_norm * max(histcounts(sprung_mass_accel_raw, 50)) / max(y_norm);
% plot(x_norm, y_norm, 'r-', 'LineWidth', 2);
% 
% % Comfort Analysis (ISO 2631)
% ax11 = nexttile(3);
% % Weighted acceleration for comfort evaluation
% weighted_accel = sprung_mass_accel_raw; % Simplified - would need proper filtering
% plot(t_sol, weighted_accel, 'g-', 'LineWidth', 1.5);
% title('Weighted Acceleration (Comfort Analysis)','FontSize',12,'FontWeight','bold');
% xlabel('Time (s)','FontSize',11); ylabel('Weighted Accel (m/s²)','FontSize',11);
% grid on; grid minor; box on;
% % Dynamic y-axis range based on data
% weighted_range = max(abs(weighted_accel));
% ylim([-weighted_range*1.1, weighted_range*1.1]);
% % Add comfort thresholds
% yline(0.315, 'r--', 'Uncomfortable', 'LineWidth', 1);
% yline(0.63, 'r-', 'Very Uncomfortable', 'LineWidth', 1);
% 
% % Simulation Summary
% ax12 = nexttile(4);
% axis off;
% % Create summary text
% summary_text = {
%     'SIMULATION SUMMARY';
%     '==================';
%     '';
%     sprintf('Vehicle Model: %s', 'Quarter Car');
%     sprintf('Scenario: %s', strrep(damage_scenario_selection,'_',' '));
%     sprintf('Road Class: %s', road_class);
%     sprintf('Speed: %.1f km/h', Vehicle_Speed);
%     sprintf('Road Length: %.0f m', Road_Length);
%     sprintf('Sampling Rate: %.0f Hz', Sampling_Freq);
%     '';
%     'PERFORMANCE METRICS';
%     '==================';
%     sprintf('Max Body Displacement: %.2f mm', max(abs(Sprung_Mass_Disp*1000)));
%     sprintf('RMS Body Displacement: %.2f mm', rms(Sprung_Mass_Disp*1000));
%     sprintf('Max Acceleration: %.2f m/s²', max(abs(sprung_mass_accel_raw)));
%     sprintf('RMS Acceleration: %.2f m/s²', rms(sprung_mass_accel_raw));
%     sprintf('Max Road Elevation: %.2f mm', max(abs(z_r_x*1000)));
%     sprintf('RMS Road Elevation: %.2f mm', rms(z_r_x*1000));
%     '';
%     'VEHICLE PARAMETERS';
%     '==================';
%     sprintf('Sprung Mass: %.0f kg', Sprung_Mass);
%     sprintf('Unsprung Mass: %.0f kg', Unsprung_Mass);
%     sprintf('Suspension Stiffness: %.0f N/m', Suspension_Stiffness);
%     sprintf('Suspension Damping: %.0f Ns/m', Suspension_Damping);
%     sprintf('Tire Stiffness: %.0f N/m', Tire_Stiffness);
%     sprintf('Tire Damping: %.0f Ns/m', Tire_Damping);
% };
% 
% text(0.05, 0.95, summary_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
%      'FontSize', 9, 'FontName', 'Courier New', 'BackgroundColor', 'white');
% 
% % Save performance metrics figure
% if save_results
%     perf_fig_path = fullfile(save_folder, 'simulation_performance_metrics.png');
%     saveas(gcf, perf_fig_path);
%     disp(['Saved: ' perf_fig_path]);
% end

%% 8.3) 3D Road Profile Visualization
% figure('Name', [figName ' - 3D Road Profile'], 'Position', [250, 200, 1200, 800]);
% 
% % Create 3D road surface
% subplot(2,2,1);
% % Generate road width for 3D effect
% road_width = 3.5; % meters (typical road width)
% y_road = linspace(-road_width/2, road_width/2, 20);
% [X_road, Y_road] = meshgrid(spatial_road, y_road);
% 
% % Create 3D road surface with elevation
% Z_road = repmat(z_r_x*1000, length(y_road), 1);
% 
% % Add some lateral variation for realism
% lateral_variation = 0.5 * sin(2*pi*spatial_road/50) .* (y_road' * ones(1,length(spatial_road)));
% Z_road = Z_road + lateral_variation;
% 
% % Plot 3D road surface
% surf(X_road, Y_road, Z_road, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
% colormap(jet);
% colorbar;
% title('3D Road Profile Surface', 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Distance (m)', 'FontSize', 11);
% ylabel('Road Width (m)', 'FontSize', 11);
% zlabel('Elevation (mm)', 'FontSize', 11);
% view(45, 30);
% grid on;
% 
% % 3D road profile with vehicle path
% subplot(2,2,2);
% % Plot road centerline
% plot3(spatial_road, zeros(size(spatial_road)), z_r_x*1000, 'b-', 'LineWidth', 3);
% hold on;
% 
% % Plot vehicle path (slightly offset for visibility)
% vehicle_path_y = 0.5; % 0.5m offset from centerline
% plot3(spatial_road, vehicle_path_y*ones(size(spatial_road)), z_r_x*1000, 'r--', 'LineWidth', 2);
% 
% % Add vehicle positions at key points
% key_points = 1:round(length(spatial_road)/10):length(spatial_road);
% for i = key_points
%     % Vehicle representation (simplified)
%     vehicle_x = spatial_road(i);
%     vehicle_y = vehicle_path_y;
%     vehicle_z = z_r_x(i)*1000;
% 
%     % Draw vehicle body (rectangle)
%     vehicle_length = 4.5; % meters
%     vehicle_width = 1.8;  % meters
%     vehicle_height = 1.5; % meters
% 
%     % Vehicle body vertices
%     vertices = [
%         vehicle_x - vehicle_length/2, vehicle_y - vehicle_width/2, vehicle_z;
%         vehicle_x + vehicle_length/2, vehicle_y - vehicle_width/2, vehicle_z;
%         vehicle_x + vehicle_length/2, vehicle_y + vehicle_width/2, vehicle_z;
%         vehicle_x - vehicle_length/2, vehicle_y + vehicle_width/2, vehicle_z;
%         vehicle_x - vehicle_length/2, vehicle_y - vehicle_width/2, vehicle_z + vehicle_height;
%         vehicle_x + vehicle_length/2, vehicle_y - vehicle_width/2, vehicle_z + vehicle_height;
%         vehicle_x + vehicle_length/2, vehicle_y + vehicle_width/2, vehicle_z + vehicle_height;
%         vehicle_x - vehicle_length/2, vehicle_y + vehicle_width/2, vehicle_z + vehicle_height
%     ];
% 
%     % Draw vehicle body
%     plot3(vertices([1:4,1],1), vertices([1:4,1],2), vertices([1:4,1],3), 'k-', 'LineWidth', 2);
%     plot3(vertices([5:8,5],1), vertices([5:8,5],2), vertices([5:8,5],3), 'k-', 'LineWidth', 2);
%     plot3([vertices(1,1), vertices(5,1)], [vertices(1,2), vertices(5,2)], [vertices(1,3), vertices(5,3)], 'k-', 'LineWidth', 2);
%     plot3([vertices(2,1), vertices(6,1)], [vertices(2,2), vertices(6,2)], [vertices(2,3), vertices(6,3)], 'k-', 'LineWidth', 2);
%     plot3([vertices(3,1), vertices(7,1)], [vertices(3,2), vertices(7,2)], [vertices(3,3), vertices(7,3)], 'k-', 'LineWidth', 2);
%     plot3([vertices(4,1), vertices(8,1)], [vertices(4,2), vertices(8,2)], [vertices(4,3), vertices(8,3)], 'k-', 'LineWidth', 2);
% end
% 
% title('3D Road Profile with Vehicle Path', 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Distance (m)', 'FontSize', 11);
% ylabel('Road Width (m)', 'FontSize', 11);
% zlabel('Elevation (mm)', 'FontSize', 11);
% legend('Road Centerline', 'Vehicle Path', 'Vehicle Positions', 'Location', 'best');
% view(45, 30);
% grid on;
% 
% % 3D road profile with damage visualization
% subplot(2,2,3);
% % Create enhanced 3D road with damage highlighting
% Z_road_enhanced = Z_road;
% 
% % Highlight damage areas (where elevation changes rapidly)
% road_gradient = abs(gradient(z_r_x*1000, spatial_road));
% damage_threshold = prctile(road_gradient, 80); % Top 20% of gradient values
% damage_mask = road_gradient > damage_threshold;
% 
% % Create damage visualization
% for i = 1:length(spatial_road)
%     if damage_mask(i)
%         Z_road_enhanced(:,i) = Z_road_enhanced(:,i) + 5; % Highlight damage areas
%     end
% end
% 
% % Plot 3D road with damage highlighting
% surf(X_road, Y_road, Z_road_enhanced, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
% colormap(jet);
% colorbar;
% title('3D Road Profile with Damage Highlighting', 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Distance (m)', 'FontSize', 11);
% ylabel('Road Width (m)', 'FontSize', 11);
% zlabel('Elevation (mm)', 'FontSize', 11);
% view(45, 30);
% grid on;
% 
% % 3D road profile with vehicle response overlay
% subplot(2,2,4);
% % Plot road profile
% plot3(spatial_road, zeros(size(spatial_road)), z_r_x*1000, 'b-', 'LineWidth', 2);
% hold on;
% 
% % Plot vehicle body displacement (scaled for visibility)
% body_disp_scaled = Sprung_Mass_Disp*1000 * 10; % Scale up for visibility
% plot3(spatial_road, ones(size(spatial_road)), body_disp_scaled, 'r-', 'LineWidth', 2);
% 
% % Plot vehicle acceleration (scaled for visibility)
% accel_scaled = sprung_mass_accel_raw * 2; % Scale for visibility
% plot3(spatial_road, 2*ones(size(spatial_road)), accel_scaled, 'g-', 'LineWidth', 2);
% 
% % Add reference lines
% plot3([spatial_road(1), spatial_road(end)], [0, 0], [0, 0], 'k--', 'LineWidth', 1);
% plot3([spatial_road(1), spatial_road(end)], [1, 1], [0, 0], 'k--', 'LineWidth', 1);
% plot3([spatial_road(1), spatial_road(end)], [2, 2], [0, 0], 'k--', 'LineWidth', 1);
% 
% title('3D Road Profile with Vehicle Response', 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Distance (m)', 'FontSize', 11);
% ylabel('Response Type', 'FontSize', 11);
% zlabel('Amplitude (scaled)', 'FontSize', 11);
% legend('Road Profile', 'Body Displacement (×10)', 'Acceleration (×2)', 'Reference Lines', 'Location', 'best');
% view(45, 30);
% grid on;
% 
% % Save 3D visualization figure
% if save_results
%     fig_3d_path = fullfile(save_folder, 'simulation_3d_road_profile.png');
%     saveas(gcf, fig_3d_path);
%     disp(['Saved: ' fig_3d_path]);
% end

%% 9) Save (structured)
if save_results
    results.meta.vehicle_model = 'quarter';
    results.meta.scenario = damage_scenario_selection;
    results.meta.road_class = road_class;
    results.meta.V_kmh = Vehicle_Speed;
    results.meta.fs = Sampling_Freq;
    results.meta.dt = dt;
    results.meta.timestamp = char(datetime("now","TimeZone","local","Format","yyyy-MM-dd_HH:mm:ss"));

    results.params = params;  % structures incl. geometry/inertias
    results.time = t_sol(:);
    results.x_m = t_sol(:) * V;
    results.spatial_road = spatial_road(:);
    results.zr.x = z_r_x(:);
    results.zr.base_x = z_r_x_base(:);   % ISO baseline (no damages)
    results.zr.t = z_r_t(:);
    results.zr.dz_dt = dz_r_t(:);
    results.zr.dt    = dt;           % time step (s)
    results.zr.dx    = V*dt;         % spatial step (m)
    results.meta.fs_spatial = 1/(V*dt);

    results.inputs.U = U;
    results.inputs.labels = wheel_labels;

    results.states.X = Y_sol;
    results.outputs.heave_m = Body_Heave_m(:);
    results.outputs.accel_clean = sprung_mass_accel_raw(:);
    results.outputs.accel_meas = sprung_mass_accel(:);

    % ISO 8608 metadata
    ISOmeta = iso8608_params(road_class, 0.1);
    results.iso8608.class       = char(ISOmeta.class);
    results.iso8608.n0_cycpm    = ISOmeta.n0;
    results.iso8608.w           = 2.0;
    results.iso8608.center_map  = ISOmeta.centers; % A..H centers (m^3)
    results.iso8608.bounds_lo   = ISOmeta.lower;   % geometric lower bounds
    results.iso8608.bounds_hi   = ISOmeta.upper;   % geometric upper bounds

    % Sensor noise metadata
    results.sensor_noise.enabled                = enable_noise;
    results.sensor_noise.noise_density_g_sqrt_hz= noise_density_g_sqrt_hz;
    results.sensor_noise.bias_drift_std_g_step  = bias_drift_std_g_per_step;
    results.sensor_noise.bias_offset_g          = bias_offset_g;

    % Save in the organized folder structure
    save(fullfile(save_folder, 'simulation_results.mat'),'results','-v7.3');
    disp(['Saved: ' fullfile(save_folder, 'simulation_results.mat')]);
end

%% 10) Helper Functions

function [spatial_road, z_r_x, time_sim, z_r_t, dz_r_t, z_r_x_base] = ...
    generateRoadProfile(Road_Length, V, dt, time_road_profile, road_class, damage_scenario)

    % Base spatial grid
    spatial_road = time_road_profile * V;
    if isempty(spatial_road)
        L = 0;
    else
        L = spatial_road(end);
    end

        % ---- ISO 8608 road roughness synthesis (frequency-domain, real elevation) ----
        % ISO model: Gd(n) = Gd(n0) * (n/n0)^(-w), with n in cycles/m, units m^3
        ISO = iso8608_params(road_class, 0.1);   % n0 = 0.1 cycles/m
        w_iso   = 2.0;                           % typical waviness exponent
        % Choose where to target at n0: class center vs band center
        use_band_center = false;                 % set true to bias to band middle
        if use_band_center
            Gd0_target = sqrt(ISO.lower_C * ISO.upper_C);
        else
            Gd0_target = ISO.center;
        end

        dx = V * dt;                             % spatial step (m)
        N  = numel(spatial_road);
        L  = max( (N-1)*dx, eps );

        if N < 8
            z_r_x = zeros(size(spatial_road));
        else
            % Frequency grid (DFT bins) in cycles/m: n_k = k/L, k = 0..N-1
            n = (0 : N - 1)' / L;
            if N >= 2
                n(1) = n(2);
            end

            % Target one-sided PSD along ISO line
            Gd = Gd0_target * (n / ISO.n0) .^ (-w_iso);
            Gd(~isfinite(Gd)) = 0;

            dn = 1/L;                           % cycles/m
            % One-sided amplitude for real signal with random phase:
            % A = sqrt(2 * Gd * dn)
            A = sqrt(2 * max(Gd,0) * dn);

            % Random phase
            phi = 2*pi*rand(N,1);
            X = A .* exp(1i*phi);

            % Enforce Hermitian symmetry → real ifft
            X(1) = 0;   % zero mean elevation
            if mod(N,2) == 0
                pos = 2:(N/2);
                X(N/2+1) = real(X(N/2+1));
                X(N:-1:N/2+2) = conj(X(2:N/2));
            else
                X(N:-1:((N+3)/2)) = conj(X(2:((N+1)/2)));
            end

            % Spatial profile (meters). Scale so Welch recovers target spectrum.
            z_r_x = real(ifft(X * N)).';
        end

        % Optional gentle high-pass to suppress ultra-long drift, keep ISO band
        if N >= 9
            fs_spatial = 1/dx; nyq = fs_spatial/2;
            n_c = 1/500;                          % ~500 m wavelength cutoff
            Wn = min(max(n_c/nyq, 1e-6), 0.99);
            [b,a] = butter(3, Wn, 'high');
            z_r_x = filtfilt(b,a, z_r_x);
        end

        % Remove tiny residual mean
        if ~isempty(z_r_x)
            z_r_x = z_r_x - mean(z_r_x);
        end

        % --- ISO-STRICT characterization & rescaling (ISO 8608) ---
        % Fractional-octave smoothing (octave -> 1/3 -> 1/12), log–log fit on ISO band,
        % then rescale so the fitted Gd(n0) hits the class centre exactly.
        if numel(z_r_x) >= 32
            n0    = 0.1;           % cycles/m (ISO reference)
            dx    = V*dt;          % spatial step
            % Target = class centre at n0 for the requested class
            Gd0_target = ISO.center;
            if use_band_center
                Gd0_target = sqrt(ISO.lower_C * ISO.upper_C);  % geometric band centre
            end
            % Welch estimate of one-sided displacement PSD Gd(n)
            nfft   = 2^nextpow2(min(N, 2^16));
            winlen = min(max(2^nextpow2(round(N/8)), 1024), N);
            win    = hanning(winlen);
            nover  = round(0.5*winlen);
            [Gd_hat, n_cycpm] = pwelch(z_r_x, win, nover, nfft, 1/dx, 'onesided');

            % ISO 8608 smoothing bands
            [nc, Gd_s] = iso8608_smooth_psd_iso(n_cycpm, Gd_hat);

            % Fit straight line log10(Gd) = a + b*log10(n) on ISO reporting band
            fit_lo = 0.011;   % cycles/m
            fit_hi = 2.828;   % cycles/m
            mfit = (nc >= fit_lo) & (nc <= fit_hi) & isfinite(Gd_s) & (Gd_s > 0);
            px = log10(nc(mfit));  py = log10(Gd_s(mfit));
            p_lin = polyfit(px, py, 1);           % y = a + b x
            b = p_lin(1);  a = p_lin(2);
            w_fit   = -b;
            Gd0_fit = 10^(a + b*log10(n0));

            % Rescale profile so the fitted Gd(n0) equals the class centre
            scale = sqrt(Gd0_target / max(Gd0_fit, eps));
            z_r_x = z_r_x * scale;

            % Verify & report measured class after rescaling (optional)
            [Gd_hat2, n2] = pwelch(z_r_x, win, nover, nfft, 1/dx, 'onesided');
            [nc2, Gd_s2]  = iso8608_smooth_psd_iso(n2, Gd_hat2);
            m2 = (nc2 >= fit_lo) & (nc2 <= fit_hi) & isfinite(Gd_s2) & (Gd_s2 > 0);
            p2 = polyfit(log10(nc2(m2)), log10(Gd_s2(m2)), 1);
            b2 = p2(1);  a2 = p2(2);
            w_fit2   = -b2;
            Gd0_meas = 10^(a2 + b2*log10(n0));
            class_meas = iso8608_class_from_value(Gd0_meas, ISO);
            fprintf('ISO 8608 (strict): w = %.2f, Gd(n0) = %.3g m^3 -> class %s\n', ...
                    w_fit2, Gd0_meas, class_meas);
        end

        z_r_x_base = z_r_x;   % ISO-8608 baseline before adding deterministic damages
    
    % ---- Damage parameter tables (Low / Med / High / Extreme) ----
    S = severityConfig();

    % ---- Scenarios ----
    name = lower(damage_scenario);
    fprintf('Generating road profile | Scenario: %s | Class %s\n', damage_scenario, road_class);

    switch name
        case 'smooth'
            % nothing added; roughness only

        case 'highway_smooth'
            % Highway conditions: very smooth with occasional minor imperfections
            % Reduce base roughness for highway conditions
            z_r_x = z_r_x * 0.3;  % 30% of normal roughness
            
            % Add occasional minor bumps (like expansion joints)
            expansion_joints = 0:50:L;  % Every 50m
            for i = 1:length(expansion_joints)
                if expansion_joints(i) > 0 && expansion_joints(i) < L
                    z_r_x = add_bump(spatial_road, z_r_x, expansion_joints(i), 0.002, 0.3);
                end
            end
            
            % Add occasional minor potholes (like utility cuts)
            utility_cuts = [0.2*L, 0.4*L, 0.6*L, 0.8*L];
            for i = 1:length(utility_cuts)
                z_r_x = add_pothole(spatial_road, z_r_x, utility_cuts(i), -0.005, 0.4);
            end

        case 'city_urban'
            % Urban driving: frequent stops, starts, and various road conditions
            % Increase base roughness for urban conditions
            z_r_x = z_r_x * 1.5;  % 150% of normal roughness
            
            % Add speed bumps
            speed_bumps = [0.15*L, 0.35*L, 0.55*L, 0.75*L, 0.95*L];
            for i = 1:length(speed_bumps)
                z_r_x = add_bump(spatial_road, z_r_x, speed_bumps(i), 0.08, 0.8);
            end
            
            % Add potholes (common in urban areas)
            potholes = [0.1*L, 0.25*L, 0.45*L, 0.65*L, 0.85*L];
            for i = 1:length(potholes)
                z_r_x = add_pothole(spatial_road, z_r_x, potholes(i), -0.04, 0.6);
            end
            
            % Add utility cuts and patches
            utility_cuts = [0.2*L, 0.4*L, 0.6*L, 0.8*L];
            for i = 1:length(utility_cuts)
                z_r_x = add_rut(spatial_road, z_r_x, utility_cuts(i), -0.008, 2.0, 0.5);
            end
            
            % Add some block cracking in older sections
            old_sections = [0.3*L, 0.7*L];
            for i = 1:length(old_sections)
                start_x = old_sections(i);
                end_x = start_x + 0.1*L;
                z_r_x = apply_periodic_cracks(spatial_road, z_r_x, start_x, end_x, 0.8, -0.003, 0.08);
            end


        case 'progressive_roughness'
            % Progressive roughness with equal 1/8 segments:
            % [0,1/8) → A, [1/8,2/8) → B, ..., [7/8,1] → H.
            % Each segment is synthesized and then ISO-strict rescaled to its class target.

            x   = spatial_road(:).';           % row vector (meters)
            N   = numel(x);
            dx  = V*dt;                         % spatial step (m)
            L   = max((N-1)*dx, eps);           % total length (m)
            n0  = ISO.n0;                       % 0.1 cyc/m (reference spatial frequency)

            % Classes A..H as a cell array of chars
            classes = cellstr(num2cell(ISO.classes.'));  % {'A','B','C','D','E','F','G','H'}

            % Equal 1/8 segments
            fbreaks = linspace(0, 1, numel(classes) + 1);  % 0,1/8,...,1

            % Initialize global profile
            z_r_x(:) = 0;

            for i = 1:numel(classes)
                xs = fbreaks(i)   * L;
                xe = fbreaks(i+1) * L;

                % Exact segment mask (no overlap/crossfade)
                mloc = (x >= xs) & (x <= xe);
                Nloc = nnz(mloc);
                if Nloc < 8, continue; end

                Lloc = max((Nloc-1)*dx, eps);
                n = (0:Nloc-1)'/Lloc; if Nloc >= 2, n(1) = n(2); end

                % ISO target for this class at n0 (centre or geometric band centre)
                ISOi = iso8608_params(classes{i}, n0);
                if exist('use_band_center','var') && use_band_center
                    Gd0_target = sqrt(ISOi.lower_C * ISOi.upper_C);
                else
                    Gd0_target = ISOi.center;
                end

                % One-sided PSD along ISO line (w_iso defined above in this function)
                Gd = Gd0_target * (n / n0).^(-w_iso);  Gd(~isfinite(Gd)) = 0;
                dn = 1 / Lloc;
                A  = sqrt(2 * max(Gd,0) * dn);

                % Random phase & Hermitian spectrum → real ifft
                phi = 2*pi*rand(Nloc,1);
                X   = A .* exp(1i*phi);  X(1) = 0;
                if mod(Nloc,2) == 0
                    X(Nloc/2+1) = real(X(Nloc/2+1));
                    X(Nloc:-1:Nloc/2+2) = conj(X(2:Nloc/2));
                else
                    X(Nloc:-1:((Nloc+3)/2)) = conj(X(2:((Nloc+1)/2)));
                end
                zloc = real(ifft(X * Nloc)).';   % row

                % Gentle high-pass (cut below ~500 m wavelength) & de-mean
                if Nloc >= 9
                    fs_spatial = 1/dx; nyq = fs_spatial/2;
                    n_c = 1/500; Wn = min(max(n_c/nyq, 1e-6), 0.99);
                    [b,a] = butter(2, Wn, 'high');
                    zloc = filtfilt(b, a, zloc);
                end
                zloc = zloc - mean(zloc);

                % ----- ISO-STRICT measure & rescale each segment -----
                % Welch PSD
                nfft   = max(256, 2^nextpow2(min(Nloc, 2^16)));
                winlen = min(max(round(Nloc/8), 128), Nloc);
                win    = hanning(winlen);
                nover  = round(0.5*winlen);
                [Gd_hat, n_cycpm] = pwelch(zloc, win, nover, nfft, 1/dx, 'onesided');

                % Fractional-octave smoothing (octave → 1/3 → 1/12)
                [nc, Gd_s] = iso8608_smooth_psd_iso(n_cycpm, Gd_hat);

                % Fit on ISO reporting band 0.011–2.828 cycles/m (robust fallback if sparse)
                fit_lo = 0.011;  fit_hi = 2.828;
                mfit = (nc >= fit_lo) & (nc <= fit_hi) & isfinite(Gd_s) & (Gd_s > 0);
                if nnz(mfit) >= 3
                    p_lin = polyfit(log10(nc(mfit)), log10(Gd_s(mfit)), 1);  % y = a + b x
                    bfit  = p_lin(1);  afit = p_lin(2);
                    Gd0_fit = 10^(afit + bfit*log10(n0));
                else
                    % Fallback: nearest Welch bin to n0
                    [~,k0] = min(abs(n_cycpm - n0));
                    Gd0_fit = max(Gd_hat(k0), eps);
                end

                % Rescale to hit target class level at n0
                zloc = zloc * sqrt(Gd0_target / max(Gd0_fit, eps));
                % ----- end ISO-STRICT per-segment -----

                % Enforce C0 continuity between segments
                idx = find(mloc);
                if i == 1
                    zloc = zloc - zloc(1);
                else
                    prev_end = z_r_x(idx(1)-1);
                    zloc = zloc + (prev_end - zloc(1));
                end

                % Place the segment
                z_r_x(mloc) = zloc;
            end

            % Final global conditioning consistent with other scenarios
            if N >= 9
                fs_spatial = 1/dx; nyq = fs_spatial/2;
                n_c = 1/500; Wn = min(max(n_c/nyq, 1e-6), 0.99);
                [b,a] = butter(2, Wn, 'high');
                z_r_x = filtfilt(b, a, z_r_x);
            end
            z_r_x = z_r_x - mean(z_r_x);

            % Baseline equals final (no deterministic damages in this scenario)
            z_r_x_base = z_r_x;

        case 'bumps_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = add_bump(spatial_road, z_r_x, centers(i), S.bump.height.(sev{i}), S.bump.width.(sev{i}));
            end

        case 'potholes_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = add_pothole(spatial_road, z_r_x, centers(i),S.pothole.depth.(sev{i}), S.pothole.width.(sev{i}));
            end

        case 'rutting_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = add_rut(spatial_road, z_r_x, centers(i), S.rut.depth.(sev{i}), S.rut.main_len.(sev{i}), ...
                    S.rut.trans.(sev{i}));
            end

        case 'crack_joints_only'  % (intentionally keeping requested name)
            centers = linspace(0.15*L, 0.85*L, 4);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = add_crack(spatial_road, z_r_x, centers(i), S.crack.depth.(sev{i}), S.crack.width.(sev{i}));
            end

        case 'block_cracking_only'
            % Split road in 4 equal bands; pitch/width/depth escalate by severity
            bands = bands4(L);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = apply_periodic_cracks(spatial_road, z_r_x, bands(i,1), bands(i,2), ...
                    S.block.pitch.(sev{i}), S.block.depth.(sev{i}), S.block.width.(sev{i}));
            end

        case 'fatigue_crocodile_only'
            bands = bands4(L);
            sev = {'low','medium','high','extreme'};
            seeds = [11 12 13 14];
            for i=1:4
                z_r_x = apply_fatigue_field(spatial_road, z_r_x, bands(i,1), bands(i,2), S.fatigue.(sev{i}), seeds(i));
            end

        case 'reflection_cracking_only'
            bands = bands4(L);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = apply_periodic_cracks(spatial_road, z_r_x, bands(i,1), bands(i,2), ...
                    S.reflect.pitch.(sev{i}), S.reflect.depth.(sev{i}), S.reflect.width.(sev{i}));
            end

        case 'corrugation_shoving_only'
            bands = bands4(L);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                s = bands(i,1) + 0.10*(bands(i,2)-bands(i,1));
                e = bands(i,2) - 0.30*(bands(i,2)-bands(i,1));
                z_r_x = add_corrugation(spatial_road, z_r_x, s, e, S.corrug.amp.(sev{i}), S.corrug.lambda.(sev{i}));
                % shoving step just after ripples
                x0 = e + 0.4;
                z_r_x = add_step(spatial_road, z_r_x, x0, S.corrug.step_h.(sev{i}), S.corrug.step_len.(sev{i}));
            end

        case 'depression_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = add_rut(spatial_road, z_r_x, centers(i), ...
                    S.depress.depth.(sev{i}), S.depress.main_len.(sev{i}), S.depress.trans.(sev{i}));
            end

        case 'swelling_only'
            centers = linspace(0.15*L, 0.85*L, 4);
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = add_rut(spatial_road, z_r_x, centers(i), ...
                    S.swell.height.(sev{i}), S.swell.main_len.(sev{i}), S.swell.trans.(sev{i}));
            end

        case 'edge_fatigue_only'
            bands = bands4(L);
            sev = {'low','medium','high','extreme'};
            seeds = [21 22 23 24];
            for i=1:4
                z_r_x = apply_edge_proxy(spatial_road, z_r_x, bands(i,1), bands(i,2), S.edge.(sev{i}), seeds(i));
            end

        case {'low_severity','medium_severity','high_severity','extreme_severity'}
            sev = parse_scenario_severity(damage_scenario);  % 'low'|'medium'|'high'|'extreme'
            z_r_x = lay_all_damages_once(spatial_road, z_r_x, 0, L, S, sev);

        case 'all_damages_mixed_severity'
            Q = [0 L/4; L/4 L/2; L/2 3*L/4; 3*L/4 L];
            sev = {'low','medium','high','extreme'};
            for i=1:4
                z_r_x = lay_all_damages_once(spatial_road, z_r_x, Q(i,1), Q(i,2), S, sev{i});
            end

        otherwise
            warning('Unknown scenario "%s" -> using Smooth.', damage_scenario);
    end

    % ---- Time-domain equivalents
    time_sim = time_road_profile;
    z_r_t = z_r_x(:);
    if numel(z_r_t)<2
        dz_r_t_raw = zeros(size(z_r_t));
    else
        dz_r_t_raw = [(z_r_t(2)-z_r_t(1))/dt; diff(z_r_t)/dt];
    end
    if numel(dz_r_t_raw)>=5
        dz_r_t = movmean(dz_r_t_raw,5);
    else
        dz_r_t = dz_r_t_raw;
    end
    dz_r_t = dz_r_t(:); z_r_t = z_r_t(:);
end

% Damage utilities & laydown

function S = severityConfig()
    % Unit helpers
    IN = 0.0254;     % inch -> meter
    FT2 = 0.092903;  % ft^2 -> m^2

    %% ----------------- BUMPS (not APCS-banded; keep original style) -----------------
    S.bump.height.low     = 0.03;
    S.bump.height.medium  = 0.05;
    S.bump.height.high    = 0.08;
    S.bump.height.extreme = 0.10;

    S.bump.width.low      = 0.60;
    S.bump.width.medium   = 1.00;
    S.bump.width.high     = 1.50;
    S.bump.width.extreme  = 1.80;

    %% ----------------- POTHOLES (APCS min size >= 6" both directions; sev by count) --
    % Generator: enforce min width and scale for realism (count handled later if you aggregate).
    min_w = 6 * IN;          % 0.1524 m
    S.pothole.width.low      = max(0.30, min_w);
    S.pothole.width.medium   = max(0.50, min_w);
    S.pothole.width.high     = max(0.80, min_w);
    S.pothole.width.extreme  = max(1.00, min_w);

    % Depths are geometric bowls (APCS doesn't specify depth bands; keep plausible):
    S.pothole.depth.low      = -0.03;
    S.pothole.depth.medium   = -0.05;
    S.pothole.depth.high     = -0.09;
    S.pothole.depth.extreme  = -0.12;

    %% ----------------- CRACKS (APCS severity by opening width; 1-D proxy here) -------
    % Your 1-D profile can't set crack opening; we emulate severity by deeper/narrow bands.
    S.crack.width.low        = 0.06;     % spatial band (m) along x, NOT opening width
    S.crack.width.medium     = 0.08;
    S.crack.width.high       = 0.10;
    S.crack.width.extreme    = 0.12;

    S.crack.depth.low        = -0.0015;  % ~1.5 mm notch
    S.crack.depth.medium     = -0.0030;  % ~3 mm
    S.crack.depth.high       = -0.0050;  % ~5 mm
    S.crack.depth.extreme    = -0.0075;  % stress test

    %% ----------------- RUTTING (APCS depth bands; meters) ----------------------------
    % Low: 0.25–0.50", Medium: 0.50–1.00", High: >1.00"
    S.rut.depth.low          = -(0.25*IN);   % -0.00635 m
    S.rut.depth.medium       = -(0.50*IN);   % -0.01270 m
    S.rut.depth.high         = -(1.00*IN);   % -0.02540 m
    S.rut.depth.extreme      = -(1.50*IN);   % -0.03810 m (beyond APCS for stress tests)

    % Shape extents (kept from your design; tune as needed)
    S.rut.main_len.low = 4.0;
    S.rut.trans.low = 1.2;

    S.rut.main_len.medium = 6.0;
    S.rut.trans.medium = 1.6;

    S.rut.main_len.high = 8.0;
    S.rut.trans.high = 2.2;

    S.rut.main_len.extreme = 10.0;
    S.rut.trans.extreme = 3.0;

    %% ----------------- BLOCK CRACKING (APCS by avg block size) -----------------------
    % Low: >9 ft^2, Medium: 1–9 ft^2, High: <1 ft^2 => map to pitch ~ block side.
    side_low_m = sqrt(9*FT2);     % ~0.914 m
    side_med_m = sqrt(4*FT2);     % ~0.609 m (mid of 1–9 ft^2)
    side_high_m = sqrt(1*FT2);    % ~0.305 m
    side_ext_m  = 0.20;           % stress case

    S.block.pitch.low = max(0.80, side_low_m);
    S.block.pitch.medium = max(0.40, side_med_m);
    S.block.pitch.high = max(0.25, side_high_m);
    S.block.pitch.extreme = side_ext_m;

    % Notch proxies for the cracks forming the blocks:
    S.block.width.low = 0.06;
    S.block.width.medium = 0.08;
    S.block.width.high = 0.10;
    S.block.width.extreme = 0.12;

    S.block.depth.low = -0.0015;
    S.block.depth.medium = -0.0030;
    S.block.depth.high = -0.0050;
    S.block.depth.extreme= -0.0075;

    %% ----------------- REFLECTION CRACKING (periodic joints) -------------------------
    % APCS classifies by crack width; we keep periodicity and reuse widths/depths above.
    S.reflect.pitch.low      = 7.0;
    S.reflect.pitch.medium   = 6.0;
    S.reflect.pitch.high     = 5.0;
    S.reflect.pitch.extreme  = 4.0;

    S.reflect.width.low      = S.block.width.low;
    S.reflect.width.medium   = S.block.width.medium;
    S.reflect.width.high     = S.block.width.high;
    S.reflect.width.extreme  = S.block.width.extreme;

    S.reflect.depth.low      = S.block.depth.low;
    S.reflect.depth.medium   = S.block.depth.medium;
    S.reflect.depth.high     = S.block.depth.high;
    S.reflect.depth.extreme  = S.block.depth.extreme;

    %% ----------------- FATIGUE / ALLIGATOR (field proxy) -----------------------------
    % Visual proxy via many small bowls; severity via count/intensity.
    S.fatigue.low    = struct('len',0.30,'n',15,'dmin',-0.006,'dmax',-0.012,'wmin',0.30,'wmax',0.60);
    S.fatigue.medium = struct('len',0.35,'n',30,'dmin',-0.008,'dmax',-0.014,'wmin',0.35,'wmax',0.70);
    S.fatigue.high   = struct('len',0.40,'n',50,'dmin',-0.010,'dmax',-0.016,'wmin',0.40,'wmax',0.80);
    S.fatigue.extreme= struct('len',0.45,'n',75,'dmin',-0.012,'dmax',-0.018,'wmin',0.45,'wmax',0.90);

    %% ----------------- CORRUGATION & SHOVING (not APCS-banded) -----------------------
    S.corrug.amp.low         = 0.004;  S.corrug.amp.medium = 0.006;
    S.corrug.amp.high        = 0.008;  S.corrug.amp.extreme= 0.010;

    S.corrug.lambda.low      = 0.80;   S.corrug.lambda.medium = 0.60;
    S.corrug.lambda.high     = 0.50;   S.corrug.lambda.extreme = 0.40;

    S.corrug.step_h.low      = 0.010;  S.corrug.step_h.medium = 0.015;
    S.corrug.step_h.high     = 0.020;  S.corrug.step_h.extreme = 0.025;

    S.corrug.step_len.low    = 0.20;   S.corrug.step_len.medium = 0.25;
    S.corrug.step_len.high   = 0.30;   S.corrug.step_len.extreme = 0.35;

    %% ----------------- LOCAL DEPRESSION / SWELLING (proxies) ------------------------
    S.depress.depth.low      = -0.010; S.depress.depth.medium = -0.015;
    S.depress.depth.high     = -0.020; S.depress.depth.extreme= -0.030;

    S.depress.main_len.low   = 2.0;    S.depress.main_len.medium = 3.0;
    S.depress.main_len.high  = 4.0;    S.depress.main_len.extreme= 5.0;

    S.depress.trans.low      = 0.40;   S.depress.trans.medium = 0.50;
    S.depress.trans.high     = 0.60;   S.depress.trans.extreme= 0.80;

    S.swell.height.low       = +0.010; S.swell.height.medium = +0.015;
    S.swell.height.high      = +0.020; S.swell.height.extreme= +0.030;

    S.swell.main_len.low     = 2.0;    S.swell.main_len.medium = 3.0;
    S.swell.main_len.high    = 4.0;    S.swell.main_len.extreme= 5.0;

    S.swell.trans.low        = 0.40;   S.swell.trans.medium = 0.50;
    S.swell.trans.high       = 0.60;   S.swell.trans.extreme= 0.80;

    %% ----------------- EDGE FATIGUE PROXY (shoulder zone) ---------------------------
    S.edge.low    = struct('len',10,'depth',-0.008,'np',4,'pdmin',-0.004,'pdmax',-0.007,'wmin',0.40,'wmax',0.60);
    S.edge.medium = struct('len',12,'depth',-0.012,'np',6,'pdmin',-0.005,'pdmax',-0.008,'wmin',0.45,'wmax',0.70);
    S.edge.high   = struct('len',15,'depth',-0.016,'np',8,'pdmin',-0.006,'pdmax',-0.010,'wmin',0.50,'wmax',0.80);
    S.edge.extreme= struct('len',18,'depth',-0.020,'np',10,'pdmin',-0.007,'pdmax',-0.012,'wmin',0.55,'wmax',0.90);
    
end

function z = lay_all_damages_once(x, z, xs, xe, S, sev)
    % Lay down 11 damage types with constant spacing in [xs,xe]
    L = xe - xs;
    if L<=0
        return
    end
    types = 11; ds = L/(types+1);

    % centers for "lumped" damages
    c = xs + (1:types)*ds;

    % 1) bump
    z = add_bump(x,z,c(1), S.bump.height.(sev), S.bump.width.(sev));

    % 2) pothole
    z = add_pothole(x,z,c(2), S.pothole.depth.(sev), S.pothole.width.(sev));

    % 3) rut (depression)
    z = add_rut(x,z,c(3), S.rut.depth.(sev), S.rut.main_len.(sev), S.rut.trans.(sev));

    % 4) a single crack/joint
    z = add_crack(x,z,c(4), S.crack.depth.(sev), S.crack.width.(sev));

    % 5) block cracking: periodic local zone in a tile
    len_tile = min(12, 0.25*L);
    xs5 = max(xs, c(5)-len_tile/2);
    xe5 = min(xe, c(5)+len_tile/2);
    z = apply_periodic_cracks(x,z,xs5,xe5, S.block.pitch.(sev), S.block.depth.(sev), S.block.width.(sev));

    % 6) fatigue/crocodile field (pothole field) in a tile
    len_tile = min(14, 0.28*L);
    xs6 = max(xs, c(6)-len_tile/2);
    xe6 = min(xe, c(6)+len_tile/2);
    z = apply_fatigue_field(x, z, xs6, xe6, S.fatigue.(sev), 7);

    % 7) reflection cracking: periodic in a tile
    len_tile = min(12, 0.25*L);
    xs7 = max(xs, c(7)-len_tile/2);
    xe7 = min(xe, c(7)+len_tile/2);
    z = apply_periodic_cracks(x,z,xs7,xe7, S.reflect.pitch.(sev), S.reflect.depth.(sev), S.reflect.width.(sev));

    % 8) corrugation & shoving
    len_tile = min(12, 0.25*L);
    xs8 = max(xs, c(8)-0.6*len_tile);
    xe8 = min(xe, c(8)+0.2*len_tile);
    z = add_corrugation(x, z, xs8, xe8, S.corrug.amp.(sev), S.corrug.lambda.(sev));
    z = add_step(x, z, xe8 + 0.4, S.corrug.step_h.(sev), S.corrug.step_len.(sev));

    % 9) local depression
    z = add_rut(x, z, c(9), S.depress.depth.(sev), S.depress.main_len.(sev), S.depress.trans.(sev));

    % 10) local swelling
    z = add_rut(x, z, c(10), S.swell.height.(sev), S.swell.main_len.(sev), S.swell.trans.(sev));

    % 11) edge proxy
    len_tile = min(15, 0.30*L);
    xsE = max(xs, c(11)-len_tile/2);
    xeE = min(xe, c(11)+len_tile/2);
    z = apply_edge_proxy(x, z, xsE, xeE, S.edge.(sev), 9);
end

function M = bands4(L)
    M = [0 L/4; L/4 L/2; L/2 3*L/4; 3*L/4 L];
end

% atomic geometry primitives
function z = add_pothole(x, z, center, depth, width)
    % Bowl (negative): full "depth" (negative) at center, 0 at edges
    if isempty(x) || width<=0, return; end
    d = abs(x - center);
    m = d <= (width/2);
    z(m) = z(m) + 0.5*depth*(1 + cos(2*pi*d(m)/width));
end

function z = add_bump(x, z, center, height, width)
    % Raised half-cosine: full +height at center, 0 at edges
    if isempty(x) || width<=0, return; end
    d = abs(x - center);
    m = d <= (width/2);
    z(m) = z(m) + 0.5*height*(1 + cos(2*pi*d(m)/width));
end

function z = add_crack(x, z, center, depth, width)
    % Narrow notch (negative): full "depth" at center
    if isempty(x) || width<=0, return; end
    d = abs(x - center);
    m = d <= (width/2);
    z(m) = z(m) + 0.5*depth*(1 + cos(2*pi*d(m)/width));
end

function z = add_rut(x, z, center, depth, main_len, trans_len)
    if isempty(x) || main_len<0 || trans_len <= 0
        return
    end
    total = main_len + 2 * trans_len;
    xs = center - total/2;
    xe = center + total/2;
    xs = max(xs, x(1));
    xe = min(xe, x(end));
    if xe<=xs
        return
    end
    % main
    if main_len>0
        ms = xs + trans_len;
        me = xe - trans_len;
        m = (x>=ms) & (x<=me);
        z(m) = z(m) + depth;
    end
    % entry
    m = (x>=xs) & (x<xs+trans_len);
    if any(m)
        alpha = (x(m)-xs)/trans_len;
        z(m) = z(m) + depth*(0.5*(1 - cos(pi*alpha)));
    end
    % exit
    m = (x>xe-trans_len) & (x<=xe);
    if any(m)
        alpha = (xe - x(m))/trans_len;
        z(m) = z(m) + depth*(0.5*(1 - cos(pi*alpha)));
    end
end

function z = add_corrugation(x, z, xs, xe, amp, lambda)
    if xe<=xs || amp==0 || lambda <= 0
        return
    end
    m = (x>=xs) & (x<=xe);
    if ~any(m)
        return
    end
    xi = x(m);
    Lseg = xe - xs;
    taper = 0.5*(1 - cos(2*pi*(xi - xs)/Lseg));  % smooth on/off
    z(m) = z(m) + amp .* sin(2*pi*(xi - xs)/lambda) .* taper;
end

function z = add_step(x, z, x0, h, trans_len)
    % Smooth ramp from x0 to x0+trans_len; +h persists after
    if trans_len<=0
        return
    end
    % In ramp
    m = (x>=x0) & (x<=x0+trans_len);
    if any(m)
        alpha = (x(m) - x0)/trans_len;
        z(m) = z(m) + h * 0.5*(1 - cos(pi*alpha));
    end
    % After ramp (hold step)
    m2 = x > (x0 + trans_len);
    z(m2) = z(m2) + h;
end

function z = apply_periodic_cracks(x,z,xs,xe,pitch,depth,width)
    if xe<=xs
        return
    end
    centers = (xs + pitch/2) : pitch : xe;
    for c = centers
        z = add_crack(x, z, c, depth, width);
    end
end

function z = apply_fatigue_field(x, z, xs, xe, P, seed)
    if xe<=xs
        return
    end
    rng(seed);
    L = xe - xs;
    n = P.n;
    zstart = xs;
    for k=1:n
        c = zstart + L * rand;
        d = P.dmin + (P.dmax - P.dmin) * rand;
        w = P.wmin + (P.wmax - P.wmin) * rand;
        z = add_pothole(x, z, c, d, 0.50 * w);  % slightly sharper inner bowls
    end
end

function z = apply_edge_proxy(x, z, xs, xe, E, seed)
    if xe<=xs
        return
    end
    rng(seed);
    center = (xs+xe)/2;
    main_len = max(0,(xe-xs)-1.0);
    z = add_rut(x, z, center, E.depth, main_len, 0.6);
    cs = linspace(xs+0.5, xe-0.5, E.np);
    for c = cs
        d = E.pdmin + (E.pdmax - E.pdmin) * rand;
        w = E.wmin + (E.wmax - E.wmin) * rand;
        z = add_pothole(x, z, c, d, w);
    end
end

% Vehicle model utilities

function p = packParams(ms_corner, mu_corner, ksc, csc, ktc, ctc)
    p.ms_corner = ms_corner;
    p.mu = mu_corner;
    p.ksc = ksc;
    p.csc = csc;
    p.ktc = ktc;
    p.ctc = ctc;
end

function [U, labels] = buildRoadInputs(spatial_road, zrx, t, V)
    % zr(t) sampled along the path; dzr/dt = (dzr/dx) * V
    x = spatial_road(:);
    zr = zrx(:);
    if numel(x) >= 2
        dzr_dt = V * gradient(zr, x);
    else
        dzr_dt = zeros(size(zr));
    end
    U = [zr dzr_dt];
    labels = {'zr','dzr'};
end

function [heave_m, ddzb] = bodyOutputs(p, Y, ~, ~)
    ys = Y(:,1);
    dys = Y(:,2);
    yus = Y(:,3);
    dyus = Y(:,4);
    ms = p.ms_corner;
    ks = p.ksc;
    cs = p.csc;
    ddzb = (-cs.*(dys - dyus) - ks.*(ys - yus)) / ms;
    heave_m = ys;
end

% Quarter model (A,B + ODE)
function [A,B] = quarter_AB(p)
    ms = p.ms_corner;
    mus = p.mu;
    ks = p.ksc;
    cs = p.csc;
    kt = p.ktc;
    ct = p.ctc;
    A = [0 1 0 0;
         -ks/ms -cs/ms  ks/ms   cs/ms;
          0 0 0 1;
          ks/mus cs/mus -(ks+kt)/mus -(cs+ct)/mus];
    B = [0 0; 0 0; 0 0; kt/mus ct/mus];
end

function dy = quarterODE(t,y,U,ts,p)
    zr = interp1(ts,U(:,1),t,'linear','extrap');
    dzr = interp1(ts,U(:,2),t,'linear','extrap');
    ys = y(1);
    dys = y(2);
    yus = y(3);
    dyus = y(4);
    ms = p.ms_corner;
    mus = p.mu;
    ks = p.ksc;
    cs = p.csc;
    kt = p.ktc;
    ct = p.ctc;
    ddys = (-cs*(dys-dyus) - ks*(ys-yus))/ms;
    ddyus = ( cs*(dys-dyus) + ks*(ys-yus) - ct*(dyus-dzr) - kt*(yus-zr) )/mus;
    dy = [dys ; ddys ; dyus ; ddyus];
end

function sev = parse_scenario_severity(scenario_name)
    % Returns 'low' | 'medium' | 'high' | 'extreme' from scenario strings like
    % 'Low_Severity', 'Medium_Severity', etc. Case/spacing robust.
    tokens = lower(split(regexprep(scenario_name,'\s+','_'),'_'));  % e.g. {"low","severity"}
    valid = {'low','medium','high','extreme'};
    hit = tokens(ismember(tokens,valid));
    if isempty(hit)
        error('Could not parse severity from scenario name: %s', scenario_name);
    end
    sev = hit{1};
end

function [VH, notes] = getVehiclePreset(name)
% Returns per-corner quarter-car parameters for a few vehicle archetypes.
% All values are per corner and in SI units.

    name = strrep(lower(name),' ','_');
    switch name
        case {'golden_car','golden','iri'}
            VH.Sprung_Mass          = 250;     % kg
            VH.Unsprung_Mass        = 35;      % kg
            VH.Suspension_Stiffness = 63.3e3;  % N/m
            VH.Suspension_Damping   = 6.53e3;  % Ns/m
            VH.Tire_Stiffness       = 653e3;   % N/m
            VH.Tire_Damping         = 0;       % Ns/m
            notes = "IRI Golden Car baseline";

        case 'sport_sedan'
            VH.Sprung_Mass          = 240;
            VH.Unsprung_Mass        = 40;
            VH.Suspension_Stiffness = 85e3;    % firmer spring
            VH.Suspension_Damping   = 9.0e3;   % higher damping
            VH.Tire_Stiffness       = 800e3;   % stiffer, low-profile tire
            VH.Tire_Damping         = 200;
            notes = "Taut body, low-profile tires";

        case 'suv'
            VH.Sprung_Mass          = 320;
            VH.Unsprung_Mass        = 45;
            VH.Suspension_Stiffness = 65e3;
            VH.Suspension_Damping   = 7.5e3;
            VH.Tire_Stiffness       = 600e3;
            VH.Tire_Damping         = 150;
            notes = "Higher mass, comfort-biased damping";

        case 'light_truck'
            VH.Sprung_Mass          = 350;
            VH.Unsprung_Mass        = 55;
            VH.Suspension_Stiffness = 95e3;    % payload capacity
            VH.Suspension_Damping   = 8.5e3;
            VH.Tire_Stiffness       = 700e3;
            VH.Tire_Damping         = 250;
            notes = "Stiffer rear, moderate damping";

        case 'ev_compact'
            VH.Sprung_Mass          = 300;     % battery mass raises sprung mass
            VH.Unsprung_Mass        = 40;
            VH.Suspension_Stiffness = 75e3;
            VH.Suspension_Damping   = 8.0e3;
            VH.Tire_Stiffness       = 720e3;
            VH.Tire_Damping         = 220;
            notes = "Heavier sprung mass, balanced ride";

        case 'luxury_sedan'
            VH.Sprung_Mass          = 330;
            VH.Unsprung_Mass        = 42;
            VH.Suspension_Stiffness = 70e3;    % softer primary ride
            VH.Suspension_Damping   = 7.0e3;   % comfort-tuned
            VH.Tire_Stiffness       = 650e3;
            VH.Tire_Damping         = 180;
            notes = "Comfort focus; softer spring & damping";

        case 'offroad_4x4'
            VH.Sprung_Mass          = 310;
            VH.Unsprung_Mass        = 60;      % robust hubs/tires
            VH.Suspension_Stiffness = 55e3;    % softer spring for travel
            VH.Suspension_Damping   = 10.0e3;  % strong damping for control
            VH.Tire_Stiffness       = 500e3;   % tall, compliant tire
            VH.Tire_Damping         = 300;
            notes = "Long travel, high damping, compliant tire";

        otherwise
            error('Unknown vehicle_preset: %s', name);
    end
end

function ISO = iso8608_params(classChar, n0)
% ISO 8608 parameters and class bands.
% - classChar: 'A'..'H'
% - n0: reference spatial frequency (cycles/m), default 0.1
% Outputs (units m^3 for PSD levels):
%   ISO.classes, ISO.centers, ISO.lower, ISO.upper
%   ISO.class, ISO.center, ISO.lower_C, ISO.upper_C, ISO.n0

    if nargin < 2, n0 = 0.1; end
    classChar = upper(classChar);
    assert(ismember(classChar,'A':'H'), 'ISO 8608 class must be A..H');

    base    = 16e-6;              % A-class center at n0=0.1 cycles/m
    centers = base * 4.^(0:7);    % A..H, each 12 dB step (×4)
    bounds  = sqrt(centers(1:end-1).*centers(2:end));

    ISO.classes = 'A':'H';
    ISO.centers = centers(:);
    ISO.lower   = [0; bounds(:)];
    ISO.upper   = [bounds(:); inf];

    idx         = double(classChar) - double('A') + 1;
    ISO.class   = classChar;
    ISO.center  = centers(idx);
    ISO.lower_C = ISO.lower(idx);
    ISO.upper_C = ISO.upper(idx);
    ISO.n0      = n0;
end

function [nc, Gs] = iso8608_smooth_psd_iso(n, G)
% ISO 8608-style smoothing of displacement PSD:
% concatenate octave, 1/3-octave, and 1/12-octave band-averaged values.
% Inputs:  n (cycles/m), G (m^3, one-sided PSD)
% Outputs: nc (band centre cycles/m), Gs (band-mean PSD, m^3)
    [nco, Go]   = local_bandavg_iso(n, G, 2.^(-9:1:-5),    1);   % octaves
    [nc3, G3]   = local_bandavg_iso(n, G, 2.^(-5:1/3:-2),  3);   % 1/3-octaves
    [nc12, G12] = local_bandavg_iso(n, G, 2.^(-2:1/12:3), 12);   % 1/12-octaves
    nc = [nco, nc3, nc12];
    Gs = [Go,  G3,  G12 ];
end

function [centres, Gav] = local_bandavg_iso(n, G, centres, B)
% Average PSD inside fractional-octave bands whose centres are "centres"
% and bandwidth is defined by B (1,3,12 → octave, third, twelfth).
    r = 2^(1/(2*B));   % half-band ratio
    Gav = nan(size(centres));
    for k = 1:numel(centres)
        lo = centres(k)/r;  hi = centres(k)*r;
        m = (n >= lo) & (n < hi);
        if any(m)
            Gav(k) = mean(G(m));   % arithmetic mean within the band
        end
    end
end

function cls = iso8608_class_from_value(Gd0, ISO)
% Classify a measured Gd(n0) against ISO A..H bands contained in ISO struct.
% Returns a single character 'A'..'H'.
    cls = ISO.classes(1);
    for i = 1:numel(ISO.classes)
        if (Gd0 >= ISO.lower(i)) && (Gd0 < ISO.upper(i))
            cls = ISO.classes(i);
            return
        end
    end
    % if above top band's upper bound, it's class H
    if Gd0 >= ISO.upper(end)
        cls = ISO.classes(end);
    end
end