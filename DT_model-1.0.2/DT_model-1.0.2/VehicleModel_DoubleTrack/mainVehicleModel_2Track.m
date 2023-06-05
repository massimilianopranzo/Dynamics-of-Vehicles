% ----------------------------------------------------------------
%% Main script for a basic simulation framework with a double track vehcile model
%  authors: 
%  rev. 1.0 Mattia Piccinini & Gastone Pietro Papini Rosati
%  rev. 2.0 Edoardo Pagot
%  date:
%  rev 1.0:    13/10/2020
%  rev 2.0:    16/05/2022
%  rev 2.1:    08/07/2022 (Biral)
%       - added Fz saturation. Correceted error in Fx
%       - initial condition is now parametric in initial speed
%       - changed the braking torque parameters to adapt to a GP2 model
% ----------------------------------------------------------------

% ----------------------------
%% Initialization
% ----------------------------
initialize_environment;

% ----------------------------
%% Load vehicle data
% ----------------------------

% test_tyre_model; % some plot to visualize the curvers resulting from the
% loaded data

vehicle_data = getVehicleDataStruct();
% pacejkaParam = loadPacejkaParam();

% ----------------------------
%% Define initial conditions for the simulation
% ----------------------------
V0 = 20/3.6; % Initial speed
X0 = loadInitialConditions(V0);

% ----------------------------
%% Simulation parameters
% ----------------------------
% 1 Fixed pedal, fixed steer
% 2 Linear pedal, fixed steer
% 3 Linear pedal, linear steer
% 4 Fixed pedal, linear steer
sim_options.sim_type = 4;
sim_options.pedal = .3;
sim_options.steer_angle = 3; % [deg]

simulationPars = getSimulationParams(); 
Ts = simulationPars.times.step_size;  % integration step for the simulation (fixed step)
T0 = simulationPars.times.t0;         % starting time of the simulation
Tf = simulationPars.times.tf;         % stop time of the simulation

% ----------------------------
%% Start Simulation
% ----------------------------
fprintf('Starting Simulation\n')
tic;
model_sim = sim('Vehicle_Model_2Track');
elapsed_time_simulation = toc;
fprintf('Simulation completed\n')
fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)

% ----------------------------
%% Post-Processing
% ----------------------------
dataAnalysis(model_sim,vehicle_data,Ts);


