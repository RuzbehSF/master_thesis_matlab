
%% 1) Simulate Quarter-Car Vehicle Dynamics Over an LTPP/ProVAL Road
% -------------------------------------------------------------------------
% PURPOSE:
%   Start of the pipeline for real-world data. Reads an elevation-vs-distance
%   CSV exported from ProVAL (LTPP data), maps it to time at a chosen speed,
%   simulates a 2-DOF quarter-car, and writes a structured MAT dataset that
%   downstream scripts already expect.
%
% CORE LOGIC:
%   1) Load CSV → detect distance/elevation columns → convert units to meters.
%   2) Resample to a uniform spatial grid with dx = V/fs, then map to time.
%   3) Build quarter-car state-space model and simulate (solver selectable).
%   4) Package results (signals, states, meta) → simulation_results.mat.
%
% INPUTS:
%   None; set parameters below (vehicle preset/overrides, CSV path, speed, fs).
%
% OUTPUTS:
%   00_Outputs/01_Simulation/SimulationData/simulation_results.mat
%   + a figure summarizing the road & responses.
% -------------------------------------------------------------------------

close all; clear; clc;

%% 0) Hyperparameters and Simulation Settings

%%% LTPP/ProVAL CSV (pick the center elevation by default)
csv_path        = "LTPP_Input/Road_Profiles/profile_9.csv";
csv_channel     = "Center Elevation_Full";   % matches your header
distance_units  = "m";                       % first column is meters
elevation_units = "cm";                      % second column is centimeters

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

% Assign required variables (names preserved)
Sprung_Mass           = VH.Sprung_Mass;          % kg
Unsprung_Mass         = VH.Unsprung_Mass;        % kg
Suspension_Stiffness  = VH.Suspension_Stiffness; % N/m
Suspension_Damping    = VH.Suspension_Damping;   % Ns/m
Tire_Stiffness        = VH.Tire_Stiffness;       % N/m
Tire_Damping          = VH.Tire_Damping;         % Ns/m

% Apply overrides if provided
fns = fieldnames(vehicle_override);
for k = 1:numel(fns)
    assignin('caller', fns{k}, vehicle_override.(fns{k}));
end

Total_Mass = 4*(Sprung_Mass + Unsprung_Mass);

%%% Simulation
Vehicle_Speed = 30;      % km/h
Sampling_Freq = 200;     % Hz (time sample rate)
noise_std_g   = 0;       % accelerometer noise stdev [g]

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

%% 1) Load CSV first, then set fs so dx matches the CSV
V  = Vehicle_Speed / 3.6;             % m/s
noise_std = noise_std_g * 9.806;      % m/s^2

% --- Read LTPP/ProVAL CSV and convert to meters ---
[x_m_raw, z_m_raw, header] = readProVALcsv(csv_path, csv_channel, distance_units, elevation_units);

% Use CSV spacing to choose a consistent sample rate unless user set Sampling_Freq
[x_m_raw, iu] = unique(x_m_raw, 'stable');  z_m_raw = z_m_raw(iu);
dx_csv = median(diff(x_m_raw(:)),'omitnan');

if isempty(Sampling_Freq) || Sampling_Freq <= 0
    Sampling_Freq = V / dx_csv;        % e.g., ~22.222/0.025 ≈ 888.9 Hz for your file
end
dt = 1 / Sampling_Freq;
dx = V * dt;

% Resample to uniform spatial grid whose dx equals V/fs
xq   = (x_m_raw(1):dx:x_m_raw(end)).';
z_r_x = interp1(x_m_raw, z_m_raw, xq, 'pchip', 'extrap');

% Canonical names used elsewhere in your pipeline
spatial_road = xq(:);
time_road_profile = (spatial_road - spatial_road(1)) / V;   % start at t=0

% Time-domain equivalents
time_sim = time_road_profile(:);
z_r_t = z_r_x(:);
if numel(z_r_t) >= 2
    dz_r_t = [(z_r_t(2)-z_r_t(1))/dt; diff(z_r_t)/dt];
else
    dz_r_t = zeros(size(z_r_t));
end
if numel(dz_r_t)>=5, dz_r_t = movmean(dz_r_t,5); end
z_r_x_base = z_r_x;   % no synthetic baseline; keep equal for compatibility

Road_Length = spatial_road(end) - spatial_road(1);

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
[~, ia] = unique(round(fsorted,2),'stable');
fu = fsorted(ia); zu = zeta(idxm(ia));
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

%% 3) Build road inputs for quarter-car model
[U, wheel_labels] = buildRoadInputs(spatial_road, z_r_x, time_sim, V);

%% 4) Build quarter-car model and simulate
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
        Ad = sysd.A; Bd = sysd.B;
        N = numel(time_sim);
        X = zeros(size(A,1), N);
        X(:,1) = X0(:);
        for k=1:N-1
            X(:,k+1) = Ad*X(:,k) + Bd*U(k,:)';
        end
        Y_sol = X.'; t_sol = time_sim(:);
    otherwise
        error('Unknown solver_method: %s', solver_method);
end

%% 5) Body outputs (heave & accel) for quarter-car
[Body_Heave_m, Body_Accel_clean] = bodyOutputs(params, Y_sol, U, t_sol);

% Keep legacy variable names for your plots/save:
Sprung_Mass_Disp      = Body_Heave_m;       % m
sprung_mass_accel_raw = Body_Accel_clean;   % m/s^2

% Map states to animation-friendly variables
Sprung_Mass_Disp   = Y_sol(:,1);  % ys (m)
Unsprung_Mass_Disp = Y_sol(:,3);  % yus (m)

%% 6) Add realistic accelerometer noise
sprung_mass_accel = sprung_mass_accel_raw + noise_std*randn(size(sprung_mass_accel_raw));

%% 7) Visualization
g0 = 9.806;
figName = 'Vehicle: QUARTER | Source: LTPP/ProVAL CSV';

save_folder = fullfile('00_Outputs', '01_Simulation', 'SimulationData');
if ~exist(save_folder, 'dir'); mkdir(save_folder); end

figure('Name', figName, 'Position', [100, 70, 980, 700]);
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

% Road (distance)
ax1 = nexttile(1);
plot(spatial_road, z_r_x*1000,'LineWidth',1);
title('Road Profile (from CSV)'); xlabel('Distance (m)'); ylabel('Elevation (mm)');
grid on; grid minor; box on; legend('Road','Location','best');
road_range = max(abs(z_r_x*1000)); span = max(road_range*1.1, 1e-6); ylim([-span, span]);

% Body heave (time)
ax2 = nexttile(2);
plot(t_sol, Sprung_Mass_Disp*1000,'LineWidth',1);
title('Body Heave (CG)'); xlabel('Time (s)'); ylabel('Displacement (mm)');
grid on; grid minor; box on;
heave_range = max(abs(Sprung_Mass_Disp*1000)); span = max(heave_range*1.1, 1e-6); ylim([-span, span]);

% Noisy accel
ax3 = nexttile(3);
yyaxis left; h1=plot(t_sol, sprung_mass_accel,'-','LineWidth',1); ylabel('m/s^2'); hold on; yline(0,'k-','LineWidth',0.5);
yyaxis right; h2=plot(t_sol, sprung_mass_accel/g0,'--','LineWidth',0.9); ylabel('g');
title('Body Vertical Acceleration (Noisy)'); xlabel('Time (s)'); grid on; grid minor; box on;
legend([h1 h2],{'Accel (m/s^2)','Accel (g)'},'Location','best');
accel_range = max(abs(sprung_mass_accel));
yyaxis left;  ylim([-max(accel_range*1.1,1e-6), max(accel_range*1.1,1e-6)]);
yyaxis right; ylim([-max(accel_range*1.1/g0,1e-9), max(accel_range*1.1/g0,1e-9)]);

% Clean accel
ax4 = nexttile(4);
plot(t_sol, sprung_mass_accel_raw,'k-','LineWidth',1); hold on; yline(0,'k-','LineWidth',0.5);
title('Body Vertical Acceleration (Clean)'); xlabel('Time (s)'); ylabel('m/s^2');
grid on; grid minor; box on;
clean_accel_range = max(abs(sprung_mass_accel_raw)); span = max(clean_accel_range*1.1, 1e-6); ylim([-span, span]);

linkaxes([ax2,ax3,ax4],'x');

if save_results
    main_fig_path = fullfile(save_folder, 'simulation_results.png');
    saveas(gcf, main_fig_path);
    disp(['Saved: ' main_fig_path]);
end

%% 8) Save (structured) — identical field layout to your generator
if save_results
    results.meta.vehicle_model = 'quarter';
    results.meta.scenario = 'LTPP_CSV';
    results.meta.road_class = 'N/A';
    results.meta.V_kmh = Vehicle_Speed;
    results.meta.fs = Sampling_Freq;
    results.meta.dt = dt;
    results.meta.timestamp = char(datetime("now","TimeZone","local","Format","yyyy-MM-dd_HH:mm:ss"));
    results.meta.source_csv = char(csv_path);
    results.meta.csv_info = header;   % distance/elevation columns + units parsed

    results.params = packParams(Sprung_Mass, Unsprung_Mass, Suspension_Stiffness, Suspension_Damping, ...
                                Tire_Stiffness, Tire_Damping);

    results.time = t_sol(:);
    results.x_m  = t_sol(:) * V;
    results.spatial_road = spatial_road(:);

    results.zr.x      = z_r_x(:);
    results.zr.base_x = z_r_x_base(:);     % keep parity with generator
    results.zr.t      = z_r_t(:);
    results.zr.dz_dt  = dz_r_t(:);
    results.zr.dt     = dt;
    results.zr.dx     = V*dt;
    results.meta.fs_spatial = 1/(V*dt);

    results.inputs.U = U;
    results.inputs.labels = wheel_labels;

    results.states.X = Y_sol;
    results.outputs.heave_m    = Body_Heave_m(:);
    results.outputs.accel_clean= sprung_mass_accel_raw(:);
    results.outputs.accel_meas = sprung_mass_accel(:);

    % ISO block kept for compatibility; content is neutral
    results.iso8608.class      = 'N/A';
    results.iso8608.n0_cycpm   = 0.1;
    results.iso8608.w          = 2.0;
    results.iso8608.center_map = [];
    results.iso8608.bounds_lo  = [];
    results.iso8608.bounds_hi  = [];

    save(fullfile(save_folder, 'simulation_results.mat'),'results','-v7.3');
    disp(['Saved: ' fullfile(save_folder, 'simulation_results.mat')]);
end

%% 9) Helper Functions (CSV I/O, vehicle model, etc.)

function [x_m, z_m, header] = readProVALcsv(csv_path, channel, dist_u, elev_u)
    % Read a ProVAL/LTPP CSV. Returns distance (m), elevation (m).
    tbl = readtable(csv_path, 'VariableNamingRule','preserve');

    % distance column
    distVar = pickVar(tbl, ["Distance","Station","Milepost","Chainage","X"]);
    if distVar == "", error('No distance-like column found in %s', csv_path); end

    % elevation column: requested channel first, else guess (but never reuse distVar)
    elevVar = "";
    if strlength(channel) > 0
        elevVar = pickVar(tbl, channel);
    end
    if elevVar == "" || elevVar == distVar
        elevVar = pickVar(tbl, ["Center Elevation_Full","Center Elevation","Elevation","Profile","Center","Z","Height"]);
    end
    if elevVar == "" || elevVar == distVar
        elevVar = firstNumericNot(tbl, distVar);
    end
    if elevVar == "" || elevVar == distVar
        error('No elevation-like column (distinct from distance) found in %s', csv_path);
    end

    x = tbl.(distVar); z = tbl.(elevVar);
    good = isfinite(x) & isfinite(z); x = x(good); z = z(good);

    % parse units from headers (e.g., 'Elevation (cm)')
    du = headerUnits(tbl, distVar); eu = headerUnits(tbl, elevVar);
    if dist_u ~= "auto", du = string(dist_u); end
    if elev_u ~= "auto", eu = string(elev_u); end

    x_m = convert_distance_to_m(x, du);
    z_m = convert_elevation_to_m(z, eu);

    header = struct('distance_var',string(distVar),'elevation_var',string(elevVar), ...
                    'distance_units',string(du),'elevation_units',string(eu));
end

function name = pickVar(tbl, candidates)
    if isstring(candidates) || ischar(candidates), candidates = string(candidates); end
    vars = string(tbl.Properties.VariableNames);
    for c = candidates
        hit = vars(strcmpi(vars,c));
        if ~isempty(hit), name = hit(1); return; end
        cont = vars(contains(lower(vars), lower(c)));
        if ~isempty(cont), name = cont(1); return; end
    end
    name = "";
end

function name = firstNumericNot(tbl, exclude)
    vars = string(tbl.Properties.VariableNames);
    for v = vars
        if v == exclude
            continue;
        end
        if isnumeric(tbl.(v))
            name = v;
            return;
        end
    end
    name = "";
end

function u = headerUnits(tbl, varName)
    txt = string(tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, varName)});
    vdesc = tbl.Properties.VariableDescriptions;
    if ~isempty(vdesc)
        idx = find(strcmp(tbl.Properties.VariableNames, varName),1);
        if ~isempty(idx) && ~isempty(vdesc{idx})
            txt = string(vdesc{idx});
        end
    end
    m = regexp(txt, "\(([^)]+)\)", "tokens", "once");
    if isempty(m)
        u = "";
        return;
    end
    token = lower(strtrim(m{1}));
    map = struct("m","m","meter","m","meters","m", ...
                 "mm","mm","millimeter","mm","millimeters","mm", ...
                 "cm","cm","centimeter","cm","centimeters","cm", ...
                 "in","in","inch","in","inches","in", ...
                 "km","km","kilometer","km","kilometers","km", ...
                 "ft","ft","feet","ft","foot","ft", ...
                 "mi","mi","mile","mi","miles","mi");
    if isfield(map, token)
        u = map.(token);
    else
        u = token;
    end
end

function x = convert_distance_to_m(x, u)
    switch lower(string(u))
        case "m"
        case "km", x = x*1000;
        case "ft", x = x*0.3048;
        case "mi", x = x*1609.344;
        otherwise % assume meters
    end
end

function z = convert_elevation_to_m(z, u)
    switch lower(string(u))
        case "m"
        case "cm", z = z/100;
        case "mm", z = z/1000;
        case "in", z = z*0.0254;
        otherwise % assume meters
    end
end

function p = packParams(ms_corner, mu_corner, ksc, csc, ktc, ctc)
    p.ms_corner = ms_corner;
    p.mu = mu_corner;
    p.ksc = ksc;
    p.csc = csc;
    p.ktc = ktc;
    p.ctc = ctc;
end

function [U, labels] = buildRoadInputs(spatial_road, zrx, ~, V)
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

function [A,B] = quarter_AB(p)
    ms = p.ms_corner; mus = p.mu; ks = p.ksc; cs = p.csc; kt = p.ktc; ct = p.ctc;
    A = [0 1 0 0;
         -ks/ms -cs/ms  ks/ms   cs/ms;
          0 0 0 1;
          ks/mus cs/mus -(ks+kt)/mus -(cs+ct)/mus];
    B = [0 0; 0 0; 0 0; kt/mus ct/mus];
end

function dy = quarterODE(t,y,U,ts,p)
    zr  = interp1(ts,U(:,1),t,'linear','extrap');
    dzr = interp1(ts,U(:,2),t,'linear','extrap');
    ys = y(1); dys = y(2); yus = y(3); dyus = y(4);
    ms = p.ms_corner; mus = p.mu; ks = p.ksc; cs = p.csc; kt = p.ktc; ct = p.ctc;
    ddys  = (-cs*(dys-dyus) - ks*(ys-yus))/ms;
    ddyus = ( cs*(dys-dyus) + ks*(ys-yus) - ct*(dyus-dzr) - kt*(yus-zr) )/mus;
    dy = [dys ; ddys ; dyus ; ddyus];
end

function [heave_m, ddzb] = bodyOutputs(p, Y, ~, ~)
    ys = Y(:,1); dys = Y(:,2); yus = Y(:,3); dyus = Y(:,4);
    ms = p.ms_corner; ks = p.ksc; cs = p.csc;
    ddzb = (-cs.*(dys - dyus) - ks.*(ys - yus)) / ms;
    heave_m = ys;
end

function [VH, notes] = getVehiclePreset(name)
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
            VH.Sprung_Mass          = 240;  VH.Unsprung_Mass = 40;
            VH.Suspension_Stiffness = 85e3; VH.Suspension_Damping = 9.0e3;
            VH.Tire_Stiffness       = 800e3; VH.Tire_Damping = 200; notes="Sport Sedan";
        case 'suv'
            VH.Sprung_Mass=320; VH.Unsprung_Mass=45; VH.Suspension_Stiffness=65e3;
            VH.Suspension_Damping=7.5e3; VH.Tire_Stiffness=600e3; VH.Tire_Damping=150; notes="SUV";
        case 'light_truck'
            VH.Sprung_Mass=350; VH.Unsprung_Mass=55; VH.Suspension_Stiffness=95e3;
            VH.Suspension_Damping=8.5e3; VH.Tire_Stiffness=700e3; VH.Tire_Damping=250; notes="Light Truck";
        case 'ev_compact'
            VH.Sprung_Mass=300; VH.Unsprung_Mass=40; VH.Suspension_Stiffness=75e3;
            VH.Suspension_Damping=8.0e3; VH.Tire_Stiffness=720e3; VH.Tire_Damping=220; notes="EV Compact";
        case 'luxury_sedan'
            VH.Sprung_Mass=330; VH.Unsprung_Mass=42; VH.Suspension_Stiffness=70e3;
            VH.Suspension_Damping=7.0e3; VH.Tire_Stiffness=650e3; VH.Tire_Damping=180; notes="Luxury Sedan";
        case 'offroad_4x4'
            VH.Sprung_Mass=310; VH.Unsprung_Mass=60; VH.Suspension_Stiffness=55e3;
            VH.Suspension_Damping=10.0e3; VH.Tire_Stiffness=500e3; VH.Tire_Damping=300; notes="Offroad 4x4";
        otherwise
            error('Unknown vehicle_preset: %s', name);
    end
end