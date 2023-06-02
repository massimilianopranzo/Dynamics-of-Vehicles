%% Perform some analysis on the data
% Tw_ij, Fx_ij, Fy_ij, Mz_ij, gamma_ij, delta_ij
extra_params = model_sim.extra_params; 

% x, y, psi, u,v, Omega, Fz_ij, delta, omega_ij, alpha_ij, kappa_ij
states = model_sim.states;
%% Load vehicle data
% ---------------------------------
Lf = vehicle_data.vehicle.Lf;  % [m] Distance between vehicle CoG and front wheels axle
Lr = vehicle_data.vehicle.Lr;  % [m] Distance between vehicle CoG and front wheels axle
L  = vehicle_data.vehicle.L;   % [m] Vehicle length
Wf = vehicle_data.vehicle.Wf;  % [m] Width of front wheels axle 
Wr = vehicle_data.vehicle.Wr;  % [m] Width of rear wheels axle                   
m  = vehicle_data.vehicle.m;   % [kg] Vehicle Mass
g  = vehicle_data.vehicle.g;   % [m/s^2] Gravitational acceleration
tau_D = vehicle_data.steering_system.tau_D;  % [-] steering system ratio (pinion-rack)
i_zz = vehicle_data.vehicle.i_zz;  % [kg*m^2] Vehicle moment of inertia around z axis

% ---------------------------------
%% Extract data from simulink model
% ---------------------------------
time_sim = states.u.time;
dt = time_sim(2)-time_sim(1);

% -----------------
% Inputs
% -----------------
ped_0      = model_sim.inputs.ped_0.data;
delta_D    = model_sim.inputs.delta_D.data;

% -----------------
% States
% -----------------
x_CoM      = states.x.data;
y_CoM      = states.y.data;
psi        = states.psi.data;
u          = states.u.data;
v          = states.v.data;
Omega      = states.Omega.data;
Fz_rr      = states.Fz_rr.data;
Fz_rl      = states.Fz_rl.data;
Fz_fr      = states.Fz_fr.data;
Fz_fl      = states.Fz_fl.data;
delta      = states.delta.data;
omega_rr   = states.omega_rr.data;
omega_rl   = states.omega_rl.data;
omega_fr   = states.omega_fr.data;
omega_fl   = states.omega_fl.data;
alpha_rr   = states.alpha_rr.data;
alpha_rl   = states.alpha_rl.data;
alpha_fr   = states.alpha_fr.data;
alpha_fl   = states.alpha_fl.data;
kappa_rr   = states.kappa_rr.data;
kappa_rl   = states.kappa_rl.data;
kappa_fr   = states.kappa_fr.data;
kappa_fl   = states.kappa_fl.data;

% -----------------
% Extra Parameters
% -----------------
Tw_rr      = extra_params.Tw_rr.data;
Tw_rl      = extra_params.Tw_rl.data;
Tw_fr      = extra_params.Tw_fr.data;
Tw_fl      = extra_params.Tw_fl.data;
Fx_rr      = extra_params.Fx_rr.data;
Fx_rl      = extra_params.Fx_rl.data;
Fx_fr      = extra_params.Fx_fr.data;
Fx_fl      = extra_params.Fx_fl.data;
Fy_rr      = extra_params.Fy_rr.data;
Fy_rl      = extra_params.Fy_rl.data;
Fy_fr      = extra_params.Fy_fr.data;
Fy_fl      = extra_params.Fy_fl.data;
Mz_rr      = extra_params.Mz_rr.data;
Mz_rl      = extra_params.Mz_rl.data;
Mz_fr      = extra_params.Mz_fr.data;
Mz_fl      = extra_params.Mz_fl.data;
gamma_rr   = extra_params.gamma_rr.data;
gamma_rl   = extra_params.gamma_rl.data;
gamma_fr   = extra_params.gamma_fr.data;
gamma_fl   = extra_params.gamma_fl.data;
delta_fr   = extra_params.delta_fr.data;
delta_fl   = extra_params.delta_fl.data;

% Chassis side slip angle beta [rad]
beta = atan(v./u);

% -----------------
% Accelerations
% -----------------
% Derivatives of u, v [m/s^2]
dot_u = diff(u)/Ts;
dot_v = diff(v)/Ts;
dot_Omega = diff(Omega)/Ts;
% Total longitudinal and lateral accelerations
Ax = dot_u(1:end) - Omega(2:end).*v(2:end);
Ay = dot_v(1:end) + Omega(2:end).*u(2:end);

% figure()
% sum_F = states.Fz_rr.Data + states.Fz_rl.Data + states.Fz_fr.Data + states.Fz_fl.Data;
% plot(states.Fz_rr.Time, sum_F)

% figure()
% w_f = 9.81*0.5*(vehicle_data.vehicle.m*vehicle_data.vehicle.Lr/(vehicle_data.vehicle.Lf+vehicle_data.vehicle.Lr));
% plot(states.Fz_fr.Time, states.Fz_fr.Data-w_f);
% hold on
% w_r = 9.81*0.5*(vehicle_data.vehicle.m*vehicle_data.vehicle.Lf/(vehicle_data.vehicle.Lf+vehicle_data.vehicle.Lr));
% plot(states.Fz_fr.Time,( states.Fz_rr.Data-w_r), '--');


% Lateral load transfer -----------------
dFzf = (extra_params.Fz_fr.Data - extra_params.Fz_fl.Data) / 2;
dFzr = (extra_params.Fz_rr.Data - extra_params.Fz_rl.Data) / 2;


% Longitudinal load transfer -----------------
dFxr = (extra_params.Fx_rr.Data - extra_params.Fx_rl.Data) / 2;
dFxf = (extra_params.Fx_fr.Data - extra_params.Fx_fl.Data) / 2;
TSA = extra_params.Mz_fr.Data + extra_params.Mz_fl.Data + extra_params.Mz_rr.Data + extra_params.Mz_rl.Data;
MZtot = dFxr * Wr + dFxf * Wf + TSA;

% Measured lateral forces -----------------
Fyf = extra_params.Fy_fr.Data + extra_params.Fy_fl.Data;
Fyr = extra_params.Fy_rr.Data + extra_params.Fy_rl.Data;

% Theoretical lateral forces -----------------
Fyf_th = m * Ay * Lr / L + i_zz * dot_Omega / L - MZtot / L;
Fyr_th = m * Ay * Lf / L - i_zz * dot_Omega / L + MZtot / L;
