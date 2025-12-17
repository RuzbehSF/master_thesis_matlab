
%% 1-a) Generate an Enhanced Quarter-Car Animation
% -------------------------------------------------------------------------
% PURPOSE:
% Creates a dynamic, interactive visualization of the quarter-car simulation
% results. It shows the vehicle model driving along the road profile, with
% moving suspension components, real-time data displays, and interactive
% playback controls.
%
% CORE LOGIC (How it works):
% 1.  Loads the 'simulation_results.mat' file to get the vehicle's time-
%     history data and road profile information.
%
% 2.  Defines all animation-specific parameters, such as the dimensions of
%     the vehicle sprite, the size of the viewing window, and video
%     recording settings.
%
% 3.  Creates the main MATLAB figure, axes for plotting, and interactive
%     UI controls (e.g., Play/Pause button, speed slider, zoom buttons).
%
% 4.  Performs an initial draw, rendering the very first frame of the
%     animation, including the road segment, the car's initial position,
%     and all text displays.
%
% 5.  Enters the main animation loop, which iterates through the time-series
%     data frame by frame. In each frame, it updates the position, size,
%     and shape of all graphics objects (car body, wheel, springs, text) to
%     create the illusion of smooth movement.
%
% 6.  Manages real-time playback speed using 'pause' commands and
%     optionally writes each rendered frame to an MP4 video file using
%     MATLAB's VideoWriter.
%
% INPUTS:
% - '01_Quarter Car Modeling/simulation_results.mat': The output file from
%   the main vehicle simulation.
%
% OUTPUTS:
% - An MP4 video file of the animation.
% - A PNG image of the final animation frame.
% - A .mat file containing performance data about the animation itself.
% - A .txt summary report.
% -------------------------------------------------------------------------

close all; clc;

% --------------- Global variable declarations ---------------
global isPlaying playbackSpeed;

% --------------- load simulation results (consistent with 01_Simulation) ---------------
sim_file = fullfile('00_Outputs','01_Simulation','SimulationData','simulation_results.mat');
assert(exist(sim_file,'file')==2, 'Missing results MAT: %s. Run the simulation script first.', sim_file);

S = load(sim_file); results = S.results;

% canonical variables (mirror names used in the animation code)
t_sol   = results.time(:);
V       = results.meta.V_kmh / 3.6;      % m/s
Vehicle_Speed = results.meta.V_kmh;      % km/h (display)
Sampling_Freq = results.meta.fs;

spatial_road  = results.spatial_road(:);
z_r_x         = results.zr.x(:);         % spatial profile
z_r_t         = results.zr.t(:);         % time-domain profile along path

Y_sol               = results.states.X;          % [ys dys yus dyus]
Sprung_Mass_Disp    = results.outputs.heave_m(:);
Unsprung_Mass_Disp  = Y_sol(:,3);

sprung_mass_accel      = results.outputs.accel_meas(:);   % noisy channel
sprung_mass_accel_raw  = results.outputs.accel_clean(:);  % clean channel

road_class                = results.meta.road_class;
damage_scenario_selection = results.meta.scenario;

if isfield(results,'meta') && isfield(results.meta,'solver_method')
    solver_method = results.meta.solver_method;
else
    solver_method = 'unknown';
end

% --------------- enhanced animation params ---------------
R_tire        = 0.22;      % [m] visual tire radius
body_len      = 1.30;      % [m] body length (quarter-car sprite)
body_h        = 0.28;      % [m] body height
body_w        = 0.45;      % [m] body width (new)
hub_w         = 0.30;      % [m] unsprung block width
hub_h         = 0.12;      % [m] unsprung block height

y_hub_base    = 0.22;                               % [m] hub ref height
y_body_base   = y_hub_base + R_tire + 0.25;         % [m] body baseline above hub

win_length_m  = 22;                                 % [m] view window
target_fps    = 30;                                 % Target frames per second
frame_stride  = max(1, round(Sampling_Freq/target_fps));    % ~30 fps
realtime      = true;                               % real-time playback
make_video    = true;                               % Enable video recording

% --------------- Create output folder ---------------
output_folder = fullfile('00_Outputs', '01_Simulation', 'Animation');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Video name with output folder
video_name = fullfile(output_folder, sprintf('quartercar_animation_%s_%s.mp4', damage_scenario_selection, road_class));
video_quality = 'high';                             % 'low', 'medium', 'high'

% --------------- enhanced data validation ---------------
assert(~isempty(t_sol) && ~isempty(Sprung_Mass_Disp), 'Loaded results appear empty.');

% Check data consistency
if length(t_sol) ~= length(Sprung_Mass_Disp)
    warning('Data length mismatch. Truncating to minimum length.');
    min_len = min(length(t_sol), length(Sprung_Mass_Disp));
    t_sol = t_sol(1:min_len);
    Sprung_Mass_Disp = Sprung_Mass_Disp(1:min_len);
    if exist('Unsprung_Mass_Disp','var')
        Unsprung_Mass_Disp = Unsprung_Mass_Disp(1:min_len);
    end
    if exist('z_r_t','var')
        z_r_t = z_r_t(1:min_len);
    end
end

% --------------- precompute ---------------
N   = min([numel(t_sol), numel(Sprung_Mass_Disp), numel(Unsprung_Mass_Disp), numel(z_r_t)]);
t_sol = t_sol(:);
ys    = Sprung_Mass_Disp(:);
yus   = Unsprung_Mass_Disp(:);
zr_t  = z_r_t(:);

% Adaptive frame skipping for performance
if N > 10000
    frame_stride = max(frame_stride, round(N/5000)); % Limit to 5000 frames max
end

% Pre-allocate arrays for better performance
N_frames = length(1:frame_stride:N);
x_car_history = zeros(N_frames, 1);
yb_history = zeros(N_frames, 1);
yh_history = zeros(N_frames, 1);
frame_times = zeros(N_frames, 1);

% y-limits (auto, with headroom)
road_min = min(z_r_x); road_max = max(z_r_x);
pad_top  = 0.35; pad_bot = 0.08;
ymin = min(road_min, min(zr_t)) - pad_bot;
ymax = max([y_body_base + max(ys) + body_h/2 + 0.18, road_max + 0.18]);

% --------------- enhanced figure + initial drawing ---------------
fig = figure('Name','Enhanced Quarter-Car Animation','Color','w', 'Position',[80 80 1400 700]);
ax = axes('Parent',fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); grid(ax,'minor');
xlabel(ax,'Distance x (m)'); ylabel(ax,'Height z (m)');
title(ax, sprintf('Enhanced Quarter-Car Animation — %s   |   V=%.1f km/h   |   fs=%g Hz   |   solver=%s', ...
    strrep(damage_scenario_selection,'_',' '), Vehicle_Speed, Sampling_Freq, solver_method));

% --------------- initialize global variables ---------------
isPlaying = true;
playbackSpeed = 1.0;

% --------------- interactive controls ---------------
% Playback controls
hPlayPause = uicontrol('Style', 'pushbutton', 'String', 'Pause', ...
    'Position', [20 20 60 30], 'Callback', @togglePlayback);
hSpeed = uicontrol('Style', 'slider', 'Min', 0.1, 'Max', 3.0, 'Value', 1.0, ...
    'Position', [90 20 100 30], 'Callback', @changeSpeed);

% Zoom controls
hZoomIn = uicontrol('Style', 'pushbutton', 'String', 'Zoom In', ...
    'Position', [200 20 60 30], 'Callback', @zoomIn);
hZoomOut = uicontrol('Style', 'pushbutton', 'String', 'Zoom Out', ...
    'Position', [270 20 60 30], 'Callback', @zoomOut);

% Control callbacks
function togglePlayback(~,~)
    global isPlaying;
    isPlaying = ~isPlaying;
    if isPlaying
        set(hPlayPause, 'String', 'Pause');
    else
        set(hPlayPause, 'String', 'Play');
    end
end

function changeSpeed(hObject,~)
    global playbackSpeed;
    playbackSpeed = get(hObject, 'Value');
end

function zoomIn(~,~)
    current_lim = get(ax, 'XLim');
    center = mean(current_lim);
    width = diff(current_lim);
    new_width = width * 0.8;
    set(ax, 'XLim', [center - new_width/2, center + new_width/2]);
end

function zoomOut(~,~)
    current_lim = get(ax, 'XLim');
    center = mean(current_lim);
    width = diff(current_lim);
    new_width = width * 1.25;
    set(ax, 'XLim', [center - new_width/2, center + new_width/2]);
end

% car at time k=1
k0 = 1;
x_car = V * t_sol(k0);
camera_offset = 5.0; % meters ahead of vehicle
xwin  = [x_car - win_length_m/2 + camera_offset, x_car + win_length_m/2 + camera_offset];

% road in window - FIXED: Proper bounds checking
idx  = (spatial_road >= xwin(1)) & (spatial_road <= xwin(2));
if ~any(idx) % edge case
    [~,i1] = min(abs(spatial_road - xwin(1)));
    [~,i2] = min(abs(spatial_road - xwin(2)));
    idx = i1:i2;
end

idx = find(idx);  % make explicit index vector (mask already bounds-safe)

% Enhanced road surface with texture - FIXED: Safe indexing
if any(idx)
    hRoad = plot(ax, spatial_road(idx), z_r_x(idx), 'k-', 'LineWidth', 2.5);
    
    % Add road markings (center line)
    center_line_y = z_r_x(idx) + 0.002; % Slightly above road surface
    hCenterLine = plot(ax, spatial_road(idx), center_line_y, 'w-', 'LineWidth', 1.0, 'LineStyle', '--');
    
    % Add road shoulders
    shoulder_width = 0.5; % meters
    hShoulder1 = plot(ax, spatial_road(idx), z_r_x(idx) - shoulder_width, 'k-', 'LineWidth', 1.0);
    hShoulder2 = plot(ax, spatial_road(idx), z_r_x(idx) + shoulder_width, 'k-', 'LineWidth', 1.0);
else
    % Fallback if no valid indices
    hRoad = plot(ax, spatial_road, z_r_x, 'k-', 'LineWidth', 2.5);
    hCenterLine = plot(ax, spatial_road, z_r_x + 0.002, 'w-', 'LineWidth', 1.0, 'LineStyle', '--');
    hShoulder1 = plot(ax, spatial_road, z_r_x - 0.5, 'k-', 'LineWidth', 1.0);
    hShoulder2 = plot(ax, spatial_road, z_r_x + 0.5, 'k-', 'LineWidth', 1.0);
end

% set limits and draw once so axes pixel size is valid
xlim(ax, xwin); ylim(ax, [ymin ymax]); drawnow;

% dynamic coords
yb = y_body_base + ys(k0);
yh = y_hub_base  + yus(k0);
zr_here = zr_t(k0);

% Enhanced vehicle body with realistic proportions
hBody = rectangle(ax, 'Position',[x_car - body_len/2, yb - body_h/2, body_len, body_h], ...
    'FaceColor',[0.2 0.6 0.9], 'EdgeColor','k', 'LineWidth',1.2);

% Add vehicle windows and details
hWindow = rectangle(ax, 'Position',[x_car - body_len*0.3, yb + body_h*0.1, body_len*0.6, body_h*0.4], ...
    'FaceColor',[0.8 0.9 1.0], 'EdgeColor','k', 'LineWidth',0.8);

% Add headlights
hHeadlight1 = rectangle(ax, 'Position',[x_car + body_len*0.4, yb + body_h*0.2, 0.08, 0.06], ...
    'FaceColor','yellow', 'EdgeColor','k', 'Curvature',[1 1]);
hHeadlight2 = rectangle(ax, 'Position',[x_car + body_len*0.4, yb - body_h*0.2, 0.08, 0.06], ...
    'FaceColor','yellow', 'EdgeColor','k', 'Curvature',[1 1]);

% hub
hHub  = rectangle(ax, 'Position',[x_car - hub_w/2,  yh - hub_h/2,  hub_w,  hub_h], ...
    'FaceColor',[0.7 0.7 0.7], 'EdgeColor','k');

% --- wheel with aspect-corrected width so it appears perfectly round ---
wheel_w = local_wheel_width(ax, xwin, ymin, ymax, R_tire);
hWheel = rectangle(ax,'Position',[x_car - wheel_w/2, yh - R_tire, wheel_w, 2*R_tire], ...
    'Curvature',[1 1], 'FaceColor',[0.2 0.2 0.2], 'EdgeColor','k');

% Enhanced suspension spring with realistic appearance
ampSusp = min(0.10, 0.20*abs((yb - body_h/2) - (yh + hub_h/2)));
[xs1,ys1] = local_spring_poly(x_car - body_len*0.22, yb - body_h/2, ...
                              x_car - body_len*0.22, yh + hub_h/2, 12, ampSusp);
hSpringSusp = plot(ax, xs1, ys1, 'Color',[0 0.2 0.8], 'LineWidth',2.0);

% Add spring coils for realism
hSpringCoils = [];
for i = 1:length(xs1)-1
    if mod(i,2) == 0
        hSpringCoils(end+1) = plot(ax, [xs1(i), xs1(i+1)], [ys1(i), ys1(i+1)], 'Color',[0 0.2 0.8], 'LineWidth',3.0);
    end
end

% Enhanced damper with realistic appearance
hDamper = plot(ax, [x_car + body_len*0.22, x_car + body_len*0.22], ...
                  [yb - body_h/2,           yh + hub_h/2], 'r-', 'LineWidth',3.0);

% Add damper details
damper_center = (yb - body_h/2 + yh + hub_h/2) / 2;
hDamperBody = rectangle(ax, 'Position',[x_car + body_len*0.22 - 0.02, damper_center - 0.05, 0.04, 0.1], ...
    'FaceColor',[0.8 0.2 0.2], 'EdgeColor','k', 'LineWidth',1.0);

% tire spring (amplitude scales with tire deflection)
ampTire = min(0.08, 0.18*abs((yh - hub_h/2) - zr_here));
[xs2,ys2] = local_spring_poly(x_car, yh - hub_h/2, x_car, zr_here, 8, ampTire);
hSpringTire = plot(ax, xs2, ys2, 'k-', 'LineWidth',1.3);

% contact point
hContact = plot(ax, x_car, zr_here, 'ro', 'MarkerFaceColor','r', 'MarkerSize',4);

% Enhanced information display
hInfo = text(ax, xwin(1)+0.5, ymin + 0.25*(ymax-ymin), ...
    sprintf('t = %.2f s   x = %.1f m   v = %.1f km/h', t_sol(k0), x_car, Vehicle_Speed), ...
    'FontSize',12, 'BackgroundColor',[1 1 1 0.8], 'Margin',3, 'EdgeColor','k');

% Add performance metrics
if exist('sprung_mass_accel','var')
    hMetrics = text(ax, xwin(2)-2.0, ymin + 0.25*(ymax-ymin), ...
        sprintf('Max Accel: %.2f m/s²\nRMS Accel: %.2f m/s²', max(abs(sprung_mass_accel)), rms(sprung_mass_accel)), ...
        'FontSize',10, 'BackgroundColor',[1 1 1 0.8], 'Margin',2, 'EdgeColor','k', 'HorizontalAlignment','right');
end

% Add progress indicator
hProgress = text(ax, xwin(1)+0.5, ymax-0.1, sprintf('Frame %d/%d (%.1f%%)', k0, N, 100*k0/N), ...
    'FontSize',10, 'BackgroundColor',[1 1 1 0.8], 'Margin',2, 'EdgeColor','k');

% Add real-time acceleration plot
accel_plot_x = xwin(1) + 0.5;
accel_plot_y = ymax - 0.3;
accel_plot_w = 3.0;
accel_plot_h = 0.2;

% Create acceleration history
accel_history = [];
if exist('sprung_mass_accel','var')
    accel_history = sprung_mass_accel(k0);
    hAccelPlot = plot(ax, accel_plot_x, accel_plot_y + accel_history*accel_plot_h/2, 'b-', 'LineWidth',1.5);
end

% Add energy bar visualization
energy_x = xwin(1) + 0.5;
energy_y = ymax - 0.5;
energy_w = 2.0;
energy_h = 0.1;

% Calculate kinetic and potential energy (use state velocities)
if exist('Sprung_Mass','var') && exist('Suspension_Stiffness','var')
    vs0 = Y_sol(k0,2);  % body velocity (dys)
    kinetic_energy   = 0.5 * Sprung_Mass * vs0^2;
    potential_energy = 0.5 * Suspension_Stiffness * (ys(k0) - yus(k0))^2;
    total_energy = kinetic_energy + potential_energy;

    % Normalize energy for visualization
    max_energy = 1000; % Joules
    energy_ratio = min(total_energy / max_energy, 1.0);

    % Draw energy bar
    hEnergyBar = rectangle(ax, 'Position', [energy_x, energy_y, energy_w*energy_ratio, energy_h], ...
        'FaceColor', 'green', 'EdgeColor', 'k', 'LineWidth', 1);
end

% --------------- enhanced video setup ---------------
if make_video
    try
        vw = VideoWriter(video_name, 'MPEG-4');
        vw.Quality = 95;  % High quality
        vw.FrameRate = target_fps;
    catch
        vw = VideoWriter(video_name, 'Motion JPEG AVI');
        vw.Quality = 90;
        vw.FrameRate = target_fps;
    end
    open(vw);
end

% --------------- enhanced animation loop ---------------
start_time = tic;
frame_idx = 1;

for k = 1:frame_stride:N
    frame_start = tic;
    
    % Check if playback is paused
    global isPlaying;
    if ~isPlaying
        pause(0.1);
        continue;
    end
    
    % car pose
    x_car = V * t_sol(k);
    yb    = y_body_base + ys(k);
    yh    = y_hub_base  + yus(k);
    zr_here = zr_t(k);

    % Dynamic camera following
    camera_offset = 5.0; % meters ahead of vehicle
    xwin = [x_car - win_length_m/2 + camera_offset, x_car + win_length_m/2 + camera_offset];
    
    % Smooth camera movement
    if k > 1
        prev_xwin = get(ax, 'XLim');
        xwin = 0.8 * prev_xwin + 0.2 * xwin; % Smooth transition
    end

    % update road within window - FIXED: Proper bounds checking
    idx  = (spatial_road >= xwin(1)) & (spatial_road <= xwin(2));
    if ~any(idx)
        [~,i1] = min(abs(spatial_road - xwin(1)));
        [~,i2] = min(abs(spatial_road - xwin(2)));
        idx = i1:i2;
    end
    
    % FIXED: Ensure idx is within bounds of both arrays
    idx = idx & (1:length(idx))' <= length(spatial_road);
    idx = idx & (1:length(idx))' <= length(z_r_x);
    
    % Only update if we have valid indices
    if any(idx)
        set(hRoad,   'XData', spatial_road(idx), 'YData', z_r_x(idx));
        
        % Update center line and shoulders
        center_line_y = z_r_x(idx) + 0.002;
        set(hCenterLine, 'XData', spatial_road(idx), 'YData', center_line_y);
        set(hShoulder1, 'XData', spatial_road(idx), 'YData', z_r_x(idx) - shoulder_width);
        set(hShoulder2, 'XData', spatial_road(idx), 'YData', z_r_x(idx) + shoulder_width);
    end
    
    set(ax, 'XLim', xwin);

    % update body, hub
    set(hBody, 'Position', [x_car - body_len/2, yb - body_h/2, body_len, body_h]);
    set(hWindow, 'Position', [x_car - body_len*0.3, yb + body_h*0.1, body_len*0.6, body_h*0.4]);
    set(hHeadlight1, 'Position', [x_car + body_len*0.4, yb + body_h*0.2, 0.08, 0.06]);
    set(hHeadlight2, 'Position', [x_car + body_len*0.4, yb - body_h*0.2, 0.08, 0.06]);
    set(hHub,  'Position', [x_car - hub_w/2,   yh - hub_h/2,   hub_w,   hub_h]);

    % update wheel width (keeps it circular under any aspect)
    wheel_w = local_wheel_width(ax, xwin, ymin, ymax, R_tire);
    set(hWheel,'Position', [x_car - wheel_w/2, yh - R_tire, wheel_w, 2*R_tire]);

    % Enhanced suspension spring update
    ampSusp = min(0.10, 0.20*abs((yb - body_h/2) - (yh + hub_h/2)));
    [xs1,ys1] = local_spring_poly(x_car - body_len*0.22, yb - body_h/2, ...
                                  x_car - body_len*0.22, yh + hub_h/2, 12, ampSusp);
    set(hSpringSusp, 'XData', xs1, 'YData', ys1);
    
    % Update spring coils
    coil_idx = 1;
    for i = 1:length(xs1)-1
        if mod(i,2) == 0 && coil_idx <= length(hSpringCoils)
            set(hSpringCoils(coil_idx), 'XData', [xs1(i), xs1(i+1)], 'YData', [ys1(i), ys1(i+1)]);
            coil_idx = coil_idx + 1;
        end
    end

    % Enhanced damper update
    set(hDamper, 'XData', [x_car + body_len*0.22, x_car + body_len*0.22], ...
                 'YData', [yb - body_h/2,         yh + hub_h/2]);
    
    % Update damper body
    damper_center = (yb - body_h/2 + yh + hub_h/2) / 2;
    set(hDamperBody, 'Position', [x_car + body_len*0.22 - 0.02, damper_center - 0.05, 0.04, 0.1]);

    % tire spring update
    ampTire = min(0.08, 0.18*abs((yh - hub_h/2) - zr_here));
    [xs2,ys2] = local_spring_poly(x_car, yh - hub_h/2, x_car, zr_here, 8, ampTire);
    set(hSpringTire,'XData', xs2, 'YData', ys2);

    % update road contact marker
    set(hContact, 'XData', x_car, 'YData', zr_here);

    % Enhanced info text update
    set(hInfo,'String',sprintf('t = %.2f s   x = %.1f m   v = %.1f km/h', t_sol(k), x_car, Vehicle_Speed), ...
              'Position',[xwin(1)+0.5, ymin + 0.25*(ymax-ymin), 0]);

    % Update progress indicator
    set(hProgress,'String',sprintf('Frame %d/%d (%.1f%%)', k, N, 100*k/N), ...
                  'Position',[xwin(1)+0.5, ymax-0.1, 0]);

    % Update acceleration plot
    if exist('sprung_mass_accel','var')
        accel_history = [accel_history, sprung_mass_accel(k)];
        if length(accel_history) > 100
            accel_history = accel_history(end-99:end);
        end
        accel_x_data = accel_plot_x + (0:length(accel_history)-1)*accel_plot_w/100;
        accel_y_data = accel_plot_y + accel_history*accel_plot_h/2;
        set(hAccelPlot, 'XData', accel_x_data, 'YData', accel_y_data);
    end

    % Update energy bar (use velocities)
    if exist('Sprung_Mass','var') && exist('Suspension_Stiffness','var')
        vsk = Y_sol(k,2);  % body velocity (dys)
        kinetic_energy   = 0.5 * Sprung_Mass * vsk^2;
        potential_energy = 0.5 * Suspension_Stiffness * (ys(k) - yus(k))^2;
        total_energy = kinetic_energy + potential_energy;
    
        max_energy = 1000; % Joules
        energy_ratio = min(total_energy / max_energy, 1.0);
    
        set(hEnergyBar, 'Position', [energy_x, energy_y, energy_w*energy_ratio, energy_h]);
    end

    % Add force vector visualization
    if k > 1 && exist('Suspension_Stiffness','var') && exist('Suspension_Damping','var')
        % Calculate forces
        spring_force = Suspension_Stiffness * (ys(k) - yus(k));
        damper_force = Suspension_Damping * ( Y_sol(k,2) - Y_sol(k,4) );  % c*(dys - dyus)
        
        % Scale forces for visualization
        force_scale = 0.001;
        
        % Spring force vector
        if exist('hSpringForce','var')
            set(hSpringForce, 'XData', x_car, 'YData', yb, 'UData', 0, 'VData', spring_force*force_scale);
        else
            hSpringForce = quiver(ax, x_car, yb, 0, spring_force*force_scale, ...
                'Color', 'blue', 'LineWidth', 2, 'MaxHeadSize', 0.3);
        end
        
        % Damper force vector
        if exist('hDamperForce','var')
            set(hDamperForce, 'XData', x_car, 'YData', yb, 'UData', 0, 'VData', damper_force*force_scale);
        else
            hDamperForce = quiver(ax, x_car, yb, 0, damper_force*force_scale, ...
                'Color', 'red', 'LineWidth', 2, 'MaxHeadSize', 0.3);
        end
    end

    % Store history for trajectory visualization
    x_car_history(frame_idx) = x_car;
    yb_history(frame_idx) = yb;
    yh_history(frame_idx) = yh;

    drawnow limitrate;

    if make_video
        writeVideo(vw, getframe(fig));
    end

    % Enhanced real-time control
    global playbackSpeed;
    if realtime
        if k > 1
            target = (t_sol(k) - t_sol(k-frame_stride)) / playbackSpeed;
            lag = target - toc;
            if lag > 0, pause(lag); end
            tic;
        else
            tic;
        end
    end
    
    % Performance monitoring
    frame_times(frame_idx) = toc(frame_start);
    frame_idx = frame_idx + 1;
    
    % Display performance info every 100 frames
    if mod(k, 100*frame_stride) == 0
        avg_frame_time = mean(frame_times(1:frame_idx-1));
        fprintf('Average frame time: %.3f ms (%.1f FPS)\n', avg_frame_time*1000, 1/avg_frame_time);
    end
end

if make_video
    close(vw);
    disp(['Saved video: ', video_name]);
end

% --------------- Save animation results ---------------
% Save final figure
final_fig_path = fullfile(output_folder, 'animation_final_frame.png');
saveas(fig, final_fig_path);
disp(['Saved final frame: ', final_fig_path]);

% Save animation data
animation_data = struct();
animation_data.meta = struct();
animation_data.meta.vehicle_model = 'quarter';
animation_data.meta.scenario = damage_scenario_selection;
animation_data.meta.road_class = road_class;
animation_data.meta.V_kmh = Vehicle_Speed;
animation_data.meta.fs = Sampling_Freq;
animation_data.meta.solver_method = solver_method;
animation_data.meta.timestamp = char(datetime("now","TimeZone","local","Format","yyyy-MM-dd_HH:mm:ss"));

animation_data.performance = struct();
animation_data.performance.total_time = toc(start_time);
animation_data.performance.avg_frame_time = mean(frame_times(1:frame_idx-1));
animation_data.performance.total_frames = frame_idx-1;
animation_data.performance.target_fps = target_fps;
animation_data.performance.actual_fps = 1/mean(frame_times(1:frame_idx-1));

animation_data.trajectory = struct();
animation_data.trajectory.x_car_history = x_car_history(1:frame_idx-1);
animation_data.trajectory.yb_history = yb_history(1:frame_idx-1);
animation_data.trajectory.yh_history = yh_history(1:frame_idx-1);
animation_data.trajectory.frame_times = frame_times(1:frame_idx-1);

% Save animation data
animation_data_path = fullfile(output_folder, 'animation_data.mat');
save(animation_data_path, 'animation_data', '-v7.3');
disp(['Saved animation data: ', animation_data_path]);

% Create summary report
summary_path = fullfile(output_folder, 'animation_summary.txt');
fid = fopen(summary_path, 'w');
fprintf(fid, 'QUARTER-CAR ANIMATION SUMMARY\n');
fprintf(fid, '============================\n\n');
fprintf(fid, 'Vehicle Model: %s\n', 'Quarter Car');
fprintf(fid, 'Scenario: %s\n', damage_scenario_selection);
fprintf(fid, 'Road Class: %s\n', road_class);
fprintf(fid, 'Speed: %.1f km/h\n', Vehicle_Speed);
fprintf(fid, 'Sampling Rate: %.0f Hz\n', Sampling_Freq);
fprintf(fid, 'Solver Method: %s\n', solver_method);
fprintf(fid, '\nANIMATION PERFORMANCE\n');
fprintf(fid, '=====================\n');
fprintf(fid, 'Total Time: %.2f seconds\n', animation_data.performance.total_time);
fprintf(fid, 'Total Frames: %d\n', animation_data.performance.total_frames);
fprintf(fid, 'Target FPS: %.1f\n', animation_data.performance.target_fps);
fprintf(fid, 'Actual FPS: %.1f\n', animation_data.performance.actual_fps);
fprintf(fid, 'Average Frame Time: %.3f ms\n', animation_data.performance.avg_frame_time*1000);
fprintf(fid, '\nOUTPUT FILES\n');
fprintf(fid, '============\n');
fprintf(fid, 'Video: %s\n', video_name);
fprintf(fid, 'Final Frame: %s\n', final_fig_path);
fprintf(fid, 'Animation Data: %s\n', animation_data_path);
fprintf(fid, 'Summary Report: %s\n', summary_path);
fclose(fid);
disp(['Saved summary report: ', summary_path]);

% Display final performance statistics
fprintf('\n========== ANIMATION COMPLETED ==========\n');
fprintf('Animation completed in %.2f seconds\n', animation_data.performance.total_time);
fprintf('Average frame time: %.3f ms (%.1f FPS)\n', animation_data.performance.avg_frame_time*1000,...
    animation_data.performance.actual_fps);
fprintf('Total frames rendered: %d\n', animation_data.performance.total_frames);
fprintf('Output folder: %s\n', output_folder);
fprintf('Video saved: %s\n', video_name);

% ---------------- helper: pixel-aware circular wheel width ----------------
function wheel_w = local_wheel_width(ax, xwin, ymin, ymax, R_tire)
    % Returns data-space width that renders as a perfect circle on screen.
    oldU = ax.Units; ax.Units = 'pixels'; pos = ax.Position; ax.Units = oldU;
    pxPerX = pos(3) / diff(xwin);
    pxPerY = pos(4) / (ymax - ymin);
    wheel_w = 2*R_tire * (pxPerY / pxPerX);
end

% ---------------- helper: make a vertical spring polyline ----------------
function [xv,yv] = local_spring_poly(x1,y1,x2,y2,nturns,amp)
% Returns a zig-zag spring between (x1,y1) top and (x2,y2) bottom.
    xm = (x1 + x2)/2; ytop = y1; ybot = y2;
    if ybot > ytop, [ytop,ybot] = deal(ybot,ytop); end
    L  = ytop - ybot;
    if L < 1e-6, xv = [xm xm]; yv = [ytop ybot]; return; end
    pad = min(0.12 * L, 0.08);
    yv = linspace(ybot+pad, ytop-pad, 2*nturns+1);
    xv = xm + amp * (-1).^(0:numel(yv)-1);
    xv = [xm, xv, xm]; yv = [ybot, yv, ytop];
end