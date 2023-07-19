% ----------------------------------------------------------------
%% Main script for a basic simulation framework with a double track vehcile model
%  authors: 
%  rev. 1.0 Mattia Piccinini & Gastone Pietro Papini Rosati
%  rev. 2.0 Edoardo Pagot
%  da
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

% pacejkaParam = loadPacejkaParam();

% ----------------------------
%% Define initial conditions for the simulation
% ----------------------------

% ----------------------------
%% Simulation parameters
% ----------------------------
sim_options.slope = 1.5; %  [deg/s]
sim_options.test_type = 2; % 1 for constant u, 2 for constant steering angle
enable_export = 0; % 1 to export the data to a .mat file
enable_plot = 1; % 1 to plot the results


simulationPars = getSimulationParams(); 
Ts = simulationPars.times.step_size;  % integration step for the simulation (fixed step)
T0 = simulationPars.times.t0;         % starting time of the simulation
Tf = simulationPars.times.tf;         % stop time of the simulation

gain = 10; % Rescale pid coefficients

if sim_options.test_type == 1
  % ----------------------------
  %% Start Simulation
  % ----------------------------
  % tau_D = 12 -> delta_D / tau_D = delta, delta_D = delta * tau_D
  % For variations NOT DELATE
  % stiffness_gain_vec = [0.8 0.9, 1, 1.1, 1.2]; % 
  % camber_vec = [ -2, -1, 1, 2 ]; % [deg]
  % toe_vec = [-2, -1, 1, 2 ]; % [deg]
  % Tf_stiffness = 40; % [s]
  % Tf_camber_vec = [ 40 40 40 40 ]; % [s]
  % Tf_toe_vec = [ 40 40 40 40 ]; % [s]

  % For single simulation 
  stiffness_gain_vec = [1]; % 
  camber_vec = [0]; % [deg]
  toe_vec = [0]; % [deg]
  Tf_stiffness = 50; % [s]
  Tf_camber_vec = [ 40 40 40 40 ]; % [s]
  Tf_toe_vec = [ 40 40 40 40 ]; % [s]
  speed_slope = 0;
  sim_options.angle = 0;
  Tinit = 5;
  V0 = 50 / 3.6; % Initial speed
  suffix = '';

  % ----------------------------
  % Stiffness variation
  % ----------------------------
  for s = 1:length(stiffness_gain_vec)
    camber = 0;
    toe = 0;
    [X0, vehicle_data]  = loadInitialConditions(V0, stiffness_gain_vec(s), camber, toe);
    fprintf('Starting Simulation\n')
    tic;
    Tf = Tf_stiffness;
    model_sim = sim('Vehicle_Model_2Track');
    elapsed_time_simulation = toc;
    fprintf('Simulation completed\n')
    fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)  
    
    name_output = strcat('s', num2str(abs(stiffness_gain_vec(s)*100)), '_t0_c0');
    dataAnalysis_script
  end

  % ----------------------------
  % Camber variation
  % ----------------------------
  if length(camber_vec) > 1
  for s = 1:length(camber_vec)
    stiffness_gain = 1;
    toe = 0;
    [X0, vehicle_data]  = loadInitialConditions(V0, stiffness_gain, camber_vec(s), toe);
    fprintf('Starting Simulation\n')
    tic;
    Tf = Tf_camber_vec(s);
    model_sim = sim('Vehicle_Model_2Track');
    elapsed_time_simulation = toc;
    fprintf('Simulation completed\n')
    fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)

    if camber_vec(s) > 0
      sign = 'p';
    else
      sign = 'n';
    end
    name_output = strcat('s100_t0_c', sign, num2str(abs(camber_vec(s))));
    dataAnalysis_script
    
  end
  end
  % ----------------------------
  % toe variation
  % ----------------------------
  if length(toe_vec) > 1
  for s = 1:length(toe_vec)
    stiffness_gain = 1;
    camber = 0;
    [X0, vehicle_data]  = loadInitialConditions(V0, stiffness_gain, camber, toe_vec(s));
    fprintf('Starting Simulation\n')
    tic;
    Tf = Tf_toe_vec(s);
    model_sim = sim('Vehicle_Model_2Track');
    elapsed_time_simulation = toc;
    fprintf('Simulation completed\n')
    fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)

    if toe_vec(s) > 0
      sign = 'p';
    else
      sign = 'n';
    end
    name_output = strcat('s100_','t', sign, num2str(abs(toe_vec(s))),'_c0');
    dataAnalysis_script

  end
  end
  save("results_constspeed.mat", "model_sim")
  % plot_saved_data;
  
elseif sim_options.test_type == 2
  % ----------------------------
  %% Start Simulation
  % ----------------------------
  % tau_D = 12 -> delta_D / tau_D = delta, delta_D = delta * tau_D
%   stiffness_gain_vec = [0.9 1 1.1]; % 
%   camber_vec = [-1 1]; % [deg]
%   toe_vec = [-1 1]; % [deg]
%   Tf_stiffness = [30 30 30]; % [s] 120
%   Tf_camber_vec = 25*ones(2,1); % [s]
%   Tf_toe_vec = 25*ones(2,1); % [s]
  stiffness_gain_vec = [1]; % 
  camber_vec = [0]; % [deg]
  toe_vec = [0]; % [deg]
  Tf_stiffness = [55]; % [s] 120
  Tf_camber_vec = 25*ones(2,1); % [s]
  Tf_toe_vec = 25*ones(2,1); % [s]
  V0 = 2 / 3.6; % Initial speed 5
  speed_slope = 0.4;
  sim_options.angle = 30; % [deg] used for fixed angle test 30
  Tinit = 5;                            % second after which start to steer
  suffix = '_conststeer';
  % 100 s simulation
  % ----------------------------
  % Stiffness variation
  % ----------------------------
  for s = 1:length(stiffness_gain_vec)
    camber = 0;
    toe = 0;
    [X0, vehicle_data]  = loadInitialConditions(V0, stiffness_gain_vec(s), camber, toe);
    fprintf('Starting Simulation\n')
    tic;
    Tf = Tf_stiffness(s);
    model_sim = sim('Vehicle_Model_2Track');
    elapsed_time_simulation = toc;
    fprintf('Simulation completed\n')
    fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)  
    
    name_output = strcat('s', num2str(abs(stiffness_gain_vec(s)*100)), '_t0_c0_conststeer');
    dataAnalysis_script
  end
  % ----------------------------
  % Camber variation
  % ----------------------------
  if length(camber_vec) > 1
  for s = 1:length(camber_vec)
    stiffness_gain = 1;
    toe = 0;
    [X0, vehicle_data]  = loadInitialConditions(V0, stiffness_gain, camber_vec(s), toe);
    fprintf('Starting Simulation\n')
    tic;
    Tf = Tf_camber_vec(s);
    model_sim = sim('Vehicle_Model_2Track');
    elapsed_time_simulation = toc;
    fprintf('Simulation completed\n')
    fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)

    if camber_vec(s) > 0
      sign = 'p';
    else
      sign = 'n';
    end
    name_output = strcat('s100_','t', sign, num2str(abs(camber_vec(s))),'_c0_conststeer');
    dataAnalysis_script

  end
  end
  % ----------------------------
  % toe variation
  % ----------------------------
  if length(toe_vec) > 1
  for s = 1:length(toe_vec)
    stiffness_gain = 1;
    camber = 0;
    [X0, vehicle_data]  = loadInitialConditions(V0, stiffness_gain, camber, toe_vec(s));
    fprintf('Starting Simulation\n')
    tic;
    Tf = Tf_toe_vec(s);
    model_sim = sim('Vehicle_Model_2Track');
    elapsed_time_simulation = toc;
    fprintf('Simulation completed\n')
    fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)

    if toe_vec(s) > 0
      sign = 'p';
    else
      sign = 'n';
    end
    name_output = strcat('s100_t0_c', sign, num2str(abs(toe_vec(s))),'_conststeer');
    dataAnalysis_script

  end
  end
save("results_consteer.mat", "model_sim")
%   plot_saved_data_conststeer;
end

% ----------------------------
%% Post-Processing
% ----------------------------
% 
% if n_sim == 1
%   dataAnalysisAxleChar(model_sim{1}, vehicle_data, Ts);
% else 
%   analysis_store = cell(1, n_sim);  
%   for i=1:n_sim
%     analysis_store{i} = dataAnalysisMulti(model_sim{i}, vehicle_data, Ts);
%   end
%   plotMultiAnalysis(analysis_store);
% end
