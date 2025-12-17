
%% ==================================================================================
%                          ROAD PROFILE RECONSTRUCTION PIPELINE
% ===================================================================================
%
% PURPOSE:
% This project simulates a vehicle's dynamic response to a generated road
% profile and then uses that response (specifically, the body's vertical
% acceleration) to reconstruct the original road profile. It includes a
% complete workflow from data generation to rigorous validation.
%
% WORKFLOW OVERVIEW:
% 1. SIMULATION: A ground-truth road is created and a quarter-car model's
%    response is simulated to generate the input data for reconstruction.
%
% 2. CORE RECONSTRUCTION: The user can select between two methods:
%    - Path A: A time-domain method using two sequential ODE inversions.
%    - Path B: A frequency-domain method using transfer function deconvolution.
%
% 3. ANALYSIS & EXPORT: The reconstructed profile is extensively validated
%    against the ground truth using quantitative metrics, visual plots, and
%    stability analysis (bootstrapping). All results are packaged.
%
% 4. ADVANCED STUDIES: The robustness of the reconstruction is tested against
%    sensor noise and inaccuracies in vehicle parameter knowledge.
%
% 5. IRI ANALYSIS: The final validation is performed using the industry-standard
%    International Roughness Index (IRI).

%% 0. Initial Setup & Cleanup
% WARNING: The 'clean_project_folders' script deletes the entire '00_Outputs' 
% folder. Comment the next line out to compare against previous results.
% Deletes the previous Outputs folder for a fresh run.
clean_outputs
clear; close all; clc;

%% 1. SIMULATION
% quarter_car_simulation:
%   - Generates ground-truth road profile and simulates vehicle response.
%   - Saves results to ./00_Outputs/01_Simulation/SimulationData/
% Optional helpers:
%   - generate_simulation_animation: Creates an MP4 animation of the simulation.
%   - perform_frequency_analysis: Plots PSDs of the simulated signals.

% === SOURCE SELECTOR FOR SECTION 1 (Simulation) ===
% "Road_Profile_Generation"
% "LTPP"
% "Real_Data"

SIM_SOURCE = "Road_Profile_Generation";

switch SIM_SOURCE
    case "Road_Profile_Generation"
        % Generates the ground-truth road and simulates the vehicle response.
        % quarter_car_simulation
        % quarter_car_simulation_SANoise
        % quarter_car_simulation_SANoise_2
        quarter_car_simulation_SANoise_3

        % Visualize Each Scenario
        % visualize_scenario

        % Validates ISO 8608 road class via PSD fit and classification.
        % check_iso8608_road_class

        % (Optional) Creates an MP4 animation of the simulation.
        % generate_animation

        % (Optional) Analyzes the frequency content of the simulated acceleration.
        % perform_frequency_analysis

    case "LTPP"
        quarter_car_LTPP;

    case "Real_Data"
        % RealData_to_MAT;
        Load_Real_Data;

    otherwise
        error("Unknown SIM_SOURCE '%s'", SIM_SOURCE);
end

%% 2. CORE RECONSTRUCTION PIPELINE

% Load simulation data and generate unified configuration struct
prepare_reconstruction_inputs

% ====== Choose Profile Reconstruction Method ======
% Methods:
%   'A': Two-Step ODE Inversion (Time Domain)
%   'B': Transfer-Function Inversion (Frequency Domain)
%   'C': Independant Component Analysis (ICA)
%   'D': Inverse State-Space + Kalman Filter and RTS Smoother
%   'E': RIVA (PID-based)
RECON_METHOD = 'D';

switch upper(RECON_METHOD)
    case 'A'
        % Path A: Time-domain ODE inversion -> Spatial mapping
        Two_Steps_ODE_Inversion
        
    case 'B'
        % Path B: Frequency-domain regularized spectral inversion (FFT)
        Transfer_Function_Inversion
        % Transfer_Function_Inversion_Optimized

    case 'C'
        % Path C: Independant Component Analysis (ICA)
        ICA
        % ICA_new
        
    case 'D'
        % Path D: Inverse State-Space + Kalman Filter (+ Adaptive) and RTS Smoother
        % InverseStateSpace_KalmanFilter_RTS
        InverseStateSpace_AdaptiveKalmanFilter_RTS
        % InverseStateSpace_AdaptiveKalmanFilter_RTS_updated

    case 'E'
        % Path E: RIVA (PID-based control approach)
        RIVA
        
    otherwise
        error('Unknown RECON_METHOD: %s (use A-E)', RECON_METHOD);
end

% profile_registration_and_metrics

%% 3. RECONSTRUCTION ANALYSIS AND EXPORT

% Assesses reconstruction stability and calculates confidence intervals via bootstrapping.
calculate_reconstruction_quality

% Consolidates all results into a final .mat package, CSV files, and a manifest.
export_reconstruction_results

% Calculates a basic set of quantitative error metrics (RMSE, MAE, etc.).
calculate_reconstruction_error_metrics

% Generates visual comparison plots (overlays, error distribution, PSDs).
generate_reconstruction_error_plots

%% 4. ADVANCED STUDIES

% Analyzes how errors in vehicle parameters affect reconstruction accuracy.
run_reconstruction_sensitivity_study

% Analyzes how sensor noise in the input acceleration affects reconstruction accuracy.
reconstruction_noise_robustness 

%% 5. IRI (INTERNATIONAL ROUGHNESS INDEX) ANALYSIS

% "Road_Profile_Generation"
% "LTPP"
% "Real_Data"

SIM_SOURCE = "Road_Profile_Generation";

switch SIM_SOURCE

    case "LTPP"
        % Import IRI that you already computed in ProVAL (true profile)
        iri_of_LTPP;
        iri_of_reconstructed_profile;
        SIM_SOURCE = "LTPP";
        iri_compare_results;

    case "Road_Profile_Generation"
        % Compute IRI from the synthetic “true” profile
        iri_of_true_profile;
        iri_of_reconstructed_profile;
        SIM_SOURCE = "Road_Profile_Generation";
        iri_compare_results;

    case "Real_Data"
        iri_of_reconstructed_profile

    otherwise
        error("Unknown SIM_SOURCE '%s' (use 'LTPP' or 'Road_Profile_Generation')", string(SIM_SOURCE));
end

%% 6. Real Data Validation

RMS_VS_IRI

PSD_Slope

validation_feature_alignment

validation_distribution_check

%% --- End of Pipeline ---
clc; close all;
fprintf('Pipeline execution complete.\n');