clear
close all
clc
load("model_sim_test2.mat")
[X0, vehicle_data] = loadInitialConditions(2/3.6, 1, 0, 0);
Tf = 50;
Ts = 1e-4;
enable_plot = 1;
enable_export = 0;
Tinit = 5;
sim_options.test_type = 2;
dataAnalysis_script