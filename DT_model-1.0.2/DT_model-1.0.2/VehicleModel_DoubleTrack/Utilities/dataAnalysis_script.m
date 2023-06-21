close all
% ----------------------------------------------------------------
%% Post-Processing and Data Analysis
% ----------------------------------------------------------------

% ---------------------------------
%% Load vehicle data
% ---------------------------------
Lf      = vehicle_data.vehicle.Lf;  % [m] Distance between vehicle CoG and front wheels axle
Lr      = vehicle_data.vehicle.Lr;  % [m] Distance between vehicle CoG and front wheels axle
L       = vehicle_data.vehicle.L;   % [m] Vehicle length
Wf      = vehicle_data.vehicle.Wf;  % [m] Width of front wheels axle 
Wr      = vehicle_data.vehicle.Wr;  % [m] Width of rear wheels axle                   
m       = vehicle_data.vehicle.m;   % [kg] Vehicle Mass
g       = vehicle_data.vehicle.g;   % [m/s^2] Gravitational acceleration
tau_D   = vehicle_data.steering_system.tau_D;  % [-] steering system ratio (pinion-rack)
h_rf    = vehicle_data.front_suspension.h_rc_f;  % [m] Height of roll center of front axle
CAx     = vehicle_data.aerodynamics.CAx;  % [-] Longitudinal stiffness of tire
CAzf    = vehicle_data.aerodynamics.CAzf;  % [-] Lateral stiffness of front tire
CAzr    = vehicle_data.aerodynamics.CAzr;  % [-] Lateral stiffness of rear tire
hg      = vehicle_data.vehicle.hGs;  % [m] Height of vehicle CoG
% ---------------------------------
%% Extract data from simulink model
% ---------------------------------
time_sim = model_sim.states.u.time;
dt = time_sim(2)-time_sim(1);

if sim_options.test_type == 2
  Tinit = 2;
end
[~, start_time] = min(abs(time_sim - Tinit)); % time when ss starts

% -----------------
% Inputs
% -----------------
ped_0      = model_sim.inputs.ped_0.data;
delta_D    = model_sim.inputs.delta_D.data;

% -----------------
% States
% -----------------
x_CoM      = model_sim.states.x.data;
y_CoM      = model_sim.states.y.data;
psi        = model_sim.states.psi.data;
u          = model_sim.states.u.data;
v          = model_sim.states.v.data;
Omega      = model_sim.states.Omega.data;
Fz_rr      = model_sim.states.Fz_rr.data;
Fz_rl      = model_sim.states.Fz_rl.data;
Fz_fr      = model_sim.states.Fz_fr.data;
Fz_fl      = model_sim.states.Fz_fl.data;
delta      = model_sim.states.delta.data;
omega_rr   = model_sim.states.omega_rr.data;
omega_rl   = model_sim.states.omega_rl.data;
omega_fr   = model_sim.states.omega_fr.data;
omega_fl   = model_sim.states.omega_fl.data;
alpha_rr   = model_sim.states.alpha_rr.data;
alpha_rl   = model_sim.states.alpha_rl.data;
alpha_fr   = model_sim.states.alpha_fr.data;
alpha_fl   = model_sim.states.alpha_fl.data;
kappa_rr   = model_sim.states.kappa_rr.data;
kappa_rl   = model_sim.states.kappa_rl.data;
kappa_fr   = model_sim.states.kappa_fr.data;
kappa_fl   = model_sim.states.kappa_fl.data;

% -----------------
% Extra Parameters
% -----------------
Tw_rr      = model_sim.extra_params.Tw_rr.data;
Tw_rl      = model_sim.extra_params.Tw_rl.data;
Tw_fr      = model_sim.extra_params.Tw_fr.data;
Tw_fl      = model_sim.extra_params.Tw_fl.data;
Fx_rr      = model_sim.extra_params.Fx_rr.data;
Fx_rl      = model_sim.extra_params.Fx_rl.data;
Fx_fr      = model_sim.extra_params.Fx_fr.data;
Fx_fl      = model_sim.extra_params.Fx_fl.data;
Fy_rr      = model_sim.extra_params.Fy_rr.data;
Fy_rl      = model_sim.extra_params.Fy_rl.data;
Fy_fr      = model_sim.extra_params.Fy_fr.data;
Fy_fl      = model_sim.extra_params.Fy_fl.data;
Mz_rr      = model_sim.extra_params.Mz_rr.data;
Mz_rl      = model_sim.extra_params.Mz_rl.data;
Mz_fr      = model_sim.extra_params.Mz_fr.data;
Mz_fl      = model_sim.extra_params.Mz_fl.data;
gamma_rr   = model_sim.extra_params.gamma_rr.data;
gamma_rl   = model_sim.extra_params.gamma_rl.data;
gamma_fr   = model_sim.extra_params.gamma_fr.data;
gamma_fl   = model_sim.extra_params.gamma_fl.data;
delta_fr   = model_sim.extra_params.delta_fr.data;
delta_fl   = model_sim.extra_params.delta_fl.data;
mu_f       = model_sim.extra_params.mu_f.data;
mu_r       = model_sim.extra_params.mu_r.data;
Fy__f       = model_sim.extra_params.Fy__f.data;
Fy__r       = model_sim.extra_params.Fy__r.data;

% Chassis side slip angle beta [rad]
beta = atan(v./u);

% -----------------
% Accelerations
% -----------------
% Derivatives of u, v [m/s^2]
dot_u = diff(u)/Ts;
dot_v = diff(v)/Ts;
% Total longitudinal and lateral accelerations
Ax = dot_u(1:end) - Omega(2:end).*v(2:end);
Ay = dot_v(1:end) + Omega(2:end).*u(2:end);
% Ax low-pass filtered signal (zero-phase digital low-pass filtering)
Wn_filter = 0.01;
[b_butt,a_butt] = butter(4,Wn_filter,'low');
Ax_filt = filtfilt(b_butt,a_butt,Ax);  
dot_u_filt = filtfilt(b_butt,a_butt,dot_u);  
% Steady state lateral acceleration
Ay_ss = Omega.*u;
% Longitudinal jerk [m/s^3]
jerk_x = diff(dot_u)/Ts;

% -----------------
% Other parameters
% -----------------
% Total CoM speed [m/s]
vG = sqrt(u.^2 + v.^2);
% Steady state and transient curvature [m]
rho_ss   = Omega./vG;
rho_tran = ((dot_v.*u(1:end-1) - dot_u.*v(1:end-1)) ./ ((vG(1:end-1)).^3)) + rho_ss(1:end-1);
% Desired sinusoidal steering angle for the equivalent single track front wheel
desired_steer_atWheel = delta_D/tau_D;
n_sim = length(model_sim);
FAxc      = CAx*u.^2; % long friction
FAz_f     = CAzf*u.^2;
FAz_r     = CAzr*u.^2;

if enable_plot
  % ---------------------------------
  %% PLOTS
  % ---------------------------------

  % ---------------------------------
  %% Plot vehicle inputs
  % ---------------------------------
  fig_input = figure('Name','Inputs','NumberTitle','off', 'Color', 'w'); clf   
  % --- pedal --- %
  ax(1) = subplot(211);
  hold on
  plot(time_sim,ped_0,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('pedal $p_0$ [-]')
  xlim([0 time_sim(end)])
  % --- delta_0 --- %
  ax(2) = subplot(212);
  plot(time_sim,delta_D,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('steering angle $\delta_D$ [deg]')
  xlim([0 time_sim(end)])
  if enable_export == 1
    export_figure(fig_input, '\fig_input.eps', 'images\');
  end

  % ---------------------------------
  %% Plot vehicle motion
  % ---------------------------------
  fig_motion = figure('Name','Motion','NumberTitle','off', 'Color', 'w'); clf   
  % --- u --- %
  ax(1) = subplot(221);
  plot(time_sim,u*3.6,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$u$ [km/h]')
  xlim([0 time_sim(end)])
  % --- v --- %
  ax(2) = subplot(222);
  plot(time_sim,v,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$v$ [m/s]')
  xlim([0 time_sim(end)])
  % --- Omega --- %
  ax(3) = subplot(223);
  plot(time_sim,Omega,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\Omega$ [rad/s]')
  xlim([0 time_sim(end)])
  % --- VG --- %
  ax(3) = subplot(224);
  plot(time_sim,sqrt(u.^2 + v.^2),'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$V_G$ [rad/s]')
  xlim([0 time_sim(end)])
  if enable_export == 1
    export_figure(fig_motion, '\fig_motion.eps', 'images\');
  end

  % ---------------------------------
  %% Plot steering angles
  % ---------------------------------
  fig_steer = figure('Name','Steer','NumberTitle','off'); clf   
  % --- delta_0 --- %
  ax(1) = subplot(221);
  plot(time_sim,delta_D,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\delta_0$ [deg]')
  xlim([0 time_sim(end)])
  % --- delta_fr --- %
  ax(2) = subplot(222);
  plot(time_sim,delta_fr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\delta_{fr}$ [deg]')
  xlim([0 time_sim(end)])
  % --- delta_fl --- %
  ax(3) = subplot(223);
  hold on
  plot(time_sim,delta_fl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\delta_{fl}$ [deg]')
  xlim([0 time_sim(end)])
  % --- comparison --- %
  ax(4) = subplot(224);
  hold on
  plot(time_sim,delta_D/tau_D,'LineWidth',2)
  plot(time_sim,delta_fr,'LineWidth',2)
  plot(time_sim,delta_fl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  legend('$\delta_D/\tau_D$','$\delta_{fr}$','$\delta_{fl}$','location','best')
  xlim([0 time_sim(end)])
  if enable_export == 1
    export_figure(fig_steer, '\fig_steer.eps', 'images\');
  end

  % -------------------------------
  %% Plot lateral tire slips and lateral forces
  % -------------------------------
  fig_lat_slip = figure('Name','Lat slips & forces','NumberTitle','off', 'Color','w'); clf
  % --- alpha_rr --- %
  ax(1) = subplot(331);
  plot(time_sim,alpha_rr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\alpha_{rr}$ [deg]')
  xlim([0 time_sim(end)])
  % --- alpha_rl --- %
  ax(2) = subplot(332);
  plot(time_sim,alpha_rl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\alpha_{rl}$ [deg]')
  xlim([0 time_sim(end)])
  % --- alpha_fr --- %
  ax(3) = subplot(333);
  plot(time_sim,alpha_fr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\alpha_{fr}$ [deg]')
  xlim([0 time_sim(end)])
  % --- alpha_fl --- %
  ax(4) = subplot(334);
  plot(time_sim,alpha_fl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\alpha_{fl}$ [deg]')
  xlim([0 time_sim(end)])
  % --- Fy_rr --- %
  ax(5) = subplot(335);
  plot(time_sim,Fy_rr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$Fy_{rr}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fy_rl --- %
  ax(6) = subplot(336);
  plot(time_sim,Fy_rl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$Fy_{rl}$ [Nm]')
  xlim([0 time_sim(end)])
  % --- Fy_fr --- %
  ax(7) = subplot(337);
  plot(time_sim,Fy_fr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$Fy_{fr}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fy_fl --- %
  ax(8) = subplot(338);
  plot(time_sim,Fy_fl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$Fy_{fl}$ [N]')
  xlim([0 time_sim(end)])
  if enable_export == 1
    export_figure(fig_lat_slip, '\fig_lat_slip.eps', 'images\');
  end

  % linkaxes(ax,'x')
  clear ax


  % ---------------------------------
  %% Plot longitudinal tire slips and longitudinal forces
  % ---------------------------------
  fig_long_slip = figure('Name','Long slips & forces','NumberTitle','off', 'Color', 'w'); clf
  % --- kappa_rr --- %
  ax(1) = subplot(331);
  plot(time_sim,kappa_rr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\kappa_{rr}$ [-]')
  xlim([0 time_sim(end)])
  % --- kappa_rl --- %
  ax(2) = subplot(332);
  plot(time_sim,kappa_rl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\kappa_{rl}$ [-]')
  xlim([0 time_sim(end)])
  % --- kappa_fr --- %
  ax(3) = subplot(333);
  plot(time_sim,kappa_fr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\kappa_{fr}$ [-]')
  xlim([0 time_sim(end)])
  % --- kappa_fl --- %
  ax(4) = subplot(334);
  plot(time_sim,kappa_fl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$\kappa_{fl}$ [-]')
  xlim([0 time_sim(end)])
  % --- Fx_rr --- %
  ax(5) = subplot(335);
  plot(time_sim,Fx_rr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$Fx_{rr}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fx_rl --- %
  ax(6) = subplot(336);
  plot(time_sim,Fx_rl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$Fx_{rl}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fx_fr --- %
  ax(7) = subplot(337);
  plot(time_sim,Fx_fr,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$Fx_{fr}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fx_fl --- %
  ax(8) = subplot(338);
  plot(time_sim,Fx_fl,'LineWidth',2)
  xline(Tinit, '--r', 'LineWidth', 2);
  grid on; box on;
  title('$Fx_{fl}$ [N]')
  xlim([0 time_sim(end)])
  if enable_export == 1
    export_figure(fig_long_slip, '\fig_long_slip.eps', 'images\');
  end
  % linkaxes(ax,'x')
  clear ax


% ---------------------------------
%% Plot wheel torques and wheel rates
% ---------------------------------
fig_wheel_rate = figure('Name','Wheel rates & torques','NumberTitle','off','Color','w'); clf
% --- omega_rr --- %
ax(1) = subplot(331);
plot(time_sim,omega_rr,'LineWidth',2)
grid on; box on;
title('$\omega_{rr}$ [rad/s]')
xlim([0 time_sim(end)])
% --- omega_rl --- %
ax(2) = subplot(332);
plot(time_sim,omega_rl,'LineWidth',2)
grid on; box on;
title('$\omega_{rl}$ [rad/s]')
xlim([0 time_sim(end)])
% --- omega_fr --- %
ax(3) = subplot(333);
plot(time_sim,omega_fr,'LineWidth',2)
grid on; box on;
title('$\omega_{fr}$ [rad/s]')
xlim([0 time_sim(end)])
% --- omega_fl --- %
ax(4) = subplot(334);
plot(time_sim,omega_fl,'LineWidth',2)
grid on; box on;
title('$\omega_{fl}$ [rad/s]')
xlim([0 time_sim(end)])
% --- Tw_rr --- %
ax(5) = subplot(335);
plot(time_sim,Tw_rr,'LineWidth',2)
grid on; box on;
title('$Tw_{rr}$ [Nm]')
xlim([0 time_sim(end)])
% --- Tw_rl --- %
ax(6) = subplot(336);
plot(time_sim,Tw_rl,'LineWidth',2)
grid on; box on;
title('$Tw_{rl}$ [Nm]')
xlim([0 time_sim(end)])
% --- Tw_fr --- %
ax(7) = subplot(337);
plot(time_sim,Tw_fr,'LineWidth',2)
grid on; box on;
title('$Tw_{fr}$ [Nm]')
xlim([0 time_sim(end)])
% --- Tw_fl --- %
ax(8) = subplot(338);
plot(time_sim,Tw_fl,'LineWidth',2)
grid on; box on;
title('$Tw_{fl}$ [Nm]')
xlim([0 time_sim(end)])
linkaxes(ax,'x')
  if enable_export == 1
    export_figure(fig_wheel_rate, '\fig_wheel_rate.eps', 'images\');
  end
clear ax

% ---------------------------------
%% Plot vertical tire loads and self-aligning torques
% ---------------------------------
fig_vert_loads = figure('Name','Vert loads & ali torques','NumberTitle','off', 'Color', 'w'); clf
% --- Fz_rr --- %
ax(1) = subplot(331);
plot(time_sim,Fz_rr,'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$Fz_{rr}$ [N]')
xlim([0 time_sim(end)])
% --- Fz_rl --- %
ax(2) = subplot(332);
plot(time_sim,Fz_rl,'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$Fz_{rl}$ [N]')
xlim([0 time_sim(end)])
% --- Fz_fr --- %
ax(3) = subplot(333);
plot(time_sim,Fz_fr,'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$Fz_{fr}$ [N]')
xlim([0 time_sim(end)])
% --- Fz_fl --- %
ax(4) = subplot(334);
plot(time_sim,Fz_fl,'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$Fz_{fl}$ [N]')
xlim([0 time_sim(end)])
% --- Mz_rr --- %
ax(5) = subplot(335);
plot(time_sim,Mz_rr,'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$Mz_{rr}$ [Nm]')
xlim([0 time_sim(end)])
% --- Mz_rl --- %
ax(6) = subplot(336);
plot(time_sim,Mz_rl,'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$Mz_{rl}$ [Nm]')
xlim([0 time_sim(end)])
% --- Mz_fr --- %
ax(7) = subplot(337);
plot(time_sim,Mz_fr,'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$Mz_{fr}$ [Nm]')
xlim([0 time_sim(end)])
% --- Mz_fl --- %
ax(8) = subplot(338);
plot(time_sim,Mz_fl,'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$Mz_{fl}$ [Nm]')
xlim([0 time_sim(end)])
if enable_export == 1
  export_figure(fig_vert_loads, '\fig_vert_loads.eps', 'images\');
end

% linkaxes(ax,'x')
clear ax


% ---------------------------------
%% Plot wheel camber
% ---------------------------------
fig_camber = figure('Name','Wheel camber','NumberTitle','off', 'Color', 'w'); clf
% --- gamma_rr --- %
ax(1) = subplot(221);
plot(time_sim,gamma_rr,'LineWidth',2)
grid on; box on;
title('$\gamma_{rr}$ [deg]')
xlim([0 time_sim(end)])
% --- gamma_rl --- %
ax(2) = subplot(222);
plot(time_sim,gamma_rl,'LineWidth',2)
grid on; box on;
title('$\gamma_{rl}$ [deg]')
xlim([0 time_sim(end)])
% --- gamma_fr --- %
ax(3) = subplot(223);
plot(time_sim,gamma_fr,'LineWidth',2)
grid on; box on;
title('$\gamma_{fr}$ [deg]')
xlim([0 time_sim(end)])
% --- gamma_fl --- %
ax(4) = subplot(224);
plot(time_sim,gamma_fl,'LineWidth',2)
grid on; box on;
title('$\gamma_{fl}$ [deg]')
xlim([0 time_sim(end)])
if enable_export == 1
  export_figure(fig_camber, '\fig_camber.eps', 'images\');
end

% linkaxes(ax,'x')
clear ax

% ---------------------------------
%% Plot accelerations, chassis side slip angle and curvature
% ---------------------------------
fig_extra = figure('Name','Extra','NumberTitle','off', 'Color', 'w'); clf
% --- ax --- %
ax(1) = subplot(221);
plot(time_sim(2:end),dot_u - Omega(2:end).*v(2:end),'LineWidth',2)
hold on
plot(time_sim(2:end),diff(u)/Ts,'--g','LineWidth',2)
plot(time_sim(2:end),Ax_filt,'-.b','LineWidth',1)
plot(time_sim(2:end),dot_u_filt,'-.r','LineWidth',1)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$a_{x}$ $[m/s^2]$')
legend('$\dot{u}-\Omega v$','$\dot{u}$','filt $\dot{u}-\Omega v$','filt $\dot{u}$','Location','southeast')
xlim([0 time_sim(end)])
% --- ay --- %
ax(2) = subplot(222);
plot(time_sim(2:end),dot_v + Omega(2:end).*u(2:end),'LineWidth',2)
hold on
plot(time_sim(2:end),Omega(2:end).*u(2:end),'--g','LineWidth',1)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$a_{y}$ $[m/s^2]$')
legend('$\dot{v}+\Omega u$','$\Omega u$','Location','best')
xlim([0 time_sim(end)])
% --- beta --- %
ax(3) = subplot(223);
plot(time_sim,rad2deg(beta),'LineWidth',2)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$\beta$ [deg]')
xlim([0 time_sim(end)])
% --- rho --- %
ax(4) = subplot(224);
plot(time_sim,rho_ss,'LineWidth',2)
hold on
plot(time_sim(1:end-1),rho_tran,'--g','LineWidth',1)
xline(Tinit, '--r', 'LineWidth', 2);
grid on; box on;
title('$\rho$ [$m^{-1}$]')
legend('$\rho_{ss}$','$\rho_{transient}$','Location','best')
xlim([0 time_sim(end)])
if enable_export == 1
  export_figure(fig_extra, '\fig_extra.eps', 'images\');
end
% linkaxes(ax,'x')
clear ax

% ---------------------------------
%% Plot vehicle pose x,y,psi
% ---------------------------------
fig_pose = figure('Name','Pose','NumberTitle','off','Color','w'); clf 
% --- x --- %
ax(1) = subplot(221);
plot(time_sim,x_CoM,'LineWidth',2)
grid on; box on;
title('$x$ [m]')
xlim([0 time_sim(end)])
% --- y --- %
ax(2) = subplot(222);
plot(time_sim,y_CoM,'LineWidth',2)
grid on; box on;
title('$y$ [m]')
xlim([0 time_sim(end)])
% --- psi --- %
ax(3) = subplot(223);
plot(time_sim,rad2deg(psi),'LineWidth',2)
grid on; box on;
title('$\psi$ [deg]')
xlim([0 time_sim(end)])

linkaxes(ax,'x')
if enable_export == 1
  export_figure(fig_pose, '\fig_pose.eps', 'images\');
end
clear ax

% -------------------------------
%% Plot G-G diagram from simulation data
% -------------------------------
fig_GGplot = figure('Name','G-G plot','NumberTitle','off','Color','w'); clf
axis equal
hold on
plot3(Ay,Ax_filt,u(1:end-1),'Color',color('purple'),'LineWidth',3)
xlabel('$a_y$ [m/s$^2$]')
ylabel('$a_x$ [m/s$^2$]')
zlabel('$u$ [m/s]')
title('G-G diagram from simulation data','FontSize',18)
grid on; box on;
if enable_export == 1
  export_figure(fig_GGplot, '\fig_GGplot.eps', 'images\');
end

% -------------------------------
%% Plot vehicle path
% -------------------------------
N = length(time_sim);
fig_real_path = figure('Name','Real Vehicle Path','NumberTitle','off','Color','w'); clf
set(gca,'fontsize',16)
hold on
axis equal
xlabel('x-coord [m]')
ylabel('y-coord [m]')
title('Real Vehicle Path','FontSize',18)
plot(x_CoM,y_CoM,'Color',color('gold'),'LineWidth',2)
for i = 1:floor(N/20):N
    rot_mat = [cos(psi(i)) -sin(psi(i)) ; sin(psi(i)) cos(psi(i))];
    pos_rr = rot_mat*[-Lr -Wr/2]';
    pos_rl = rot_mat*[-Lr +Wr/2]';
    pos_fr = rot_mat*[+Lf -Wf/2]';
    pos_fl = rot_mat*[+Lf +Wf/2]';
    pos = [pos_rr pos_rl pos_fl pos_fr];
    p = patch(x_CoM(i) + pos(1,:),y_CoM(i) + pos(2,:),'blue');
    quiver(x_CoM(i), y_CoM(i), u(i)*cos(psi(i)), u(i)*sin(psi(i)), 'color', [1,0,0]);
    quiver(x_CoM(i), y_CoM(i), -v(i)*sin(psi(i)), v(i)*cos(psi(i)), 'color', [0.23,0.37,0.17]);
end
grid on; box on;
hold off
if enable_export == 1
  export_figure(fig_real_path, '\fig_real_path.eps', 'images\');
end


% % ---
% %% FORCE AS FUNCTION OF THE SLIP
% % --
%  % -------------------------------
% %% Plot lateral tire slips and lateral forces
% % -------------------------------
% figure('Name','Force slip','NumberTitle','off'), clf
% % --- alpha_rr --- %
% ax(1) = subplot(224);
% plot(alpha_rr,Fy_rr,'LineWidth',2)
% grid on; box on;
% title('$F_{rr}$ [deg]')
% % --- alpha_rl --- %
% ax(2) = subplot(223);
% plot(alpha_rl,Fy_rl,'LineWidth',2)
% grid on; box on;
% title('$F_{rl}$ [deg]')
% % --- alpha_fr --- %
% ax(3) = subplot(222);
% plot(alpha_fr,Fy_fr,'LineWidth',2)
% grid on; box on;
% title('$F_{fr}$ [deg]')
% % --- alpha_fl --- %
% ax(4) = subplot(221);
% plot(alpha_fl,Fy_fl,'LineWidth',2)
% grid on; box on;
% title('$F_{fl}$ [deg]')
% sgtitle('Lateral load $F_{y}$')
% % linkaxes(ax,'x')
% clear ax

end



%% REMOVE THE TRANSIENT
ped_0      = model_sim.inputs.ped_0.data(start_time:end);
delta_D    = model_sim.inputs.delta_D.data(start_time:end);
x_CoM      = model_sim.states.x.data(start_time:end);
y_CoM      = model_sim.states.y.data(start_time:end);
psi        = model_sim.states.psi.data(start_time:end);
u          = model_sim.states.u.data(start_time:end);
v          = model_sim.states.v.data(start_time:end);
Omega      = model_sim.states.Omega.data(start_time:end);
Fz_rr      = model_sim.states.Fz_rr.data(start_time:end);
Fz_rl      = model_sim.states.Fz_rl.data(start_time:end);
Fz_fr      = model_sim.states.Fz_fr.data(start_time:end);
Fz_fl      = model_sim.states.Fz_fl.data(start_time:end);
delta      = model_sim.states.delta.data(start_time:end);
omega_rr   = model_sim.states.omega_rr.data(start_time:end);
omega_rl   = model_sim.states.omega_rl.data(start_time:end);
omega_fr   = model_sim.states.omega_fr.data(start_time:end);
omega_fl   = model_sim.states.omega_fl.data(start_time:end);
alpha_rr   = model_sim.states.alpha_rr.data(start_time:end);
alpha_rl   = model_sim.states.alpha_rl.data(start_time:end);
alpha_fr   = model_sim.states.alpha_fr.data(start_time:end);
alpha_fl   = model_sim.states.alpha_fl.data(start_time:end);
kappa_rr   = model_sim.states.kappa_rr.data(start_time:end);
kappa_rl   = model_sim.states.kappa_rl.data(start_time:end);
kappa_fr   = model_sim.states.kappa_fr.data(start_time:end);
kappa_fl   = model_sim.states.kappa_fl.data(start_time:end);
Tw_rr      = model_sim.extra_params.Tw_rr.data(start_time:end);
Tw_rl      = model_sim.extra_params.Tw_rl.data(start_time:end);
Tw_fr      = model_sim.extra_params.Tw_fr.data(start_time:end);
Tw_fl      = model_sim.extra_params.Tw_fl.data(start_time:end);
Fx_rr      = model_sim.extra_params.Fx_rr.data(start_time:end);
Fx_rl      = model_sim.extra_params.Fx_rl.data(start_time:end);
Fx_fr      = model_sim.extra_params.Fx_fr.data(start_time:end);
Fx_fl      = model_sim.extra_params.Fx_fl.data(start_time:end);
Fy_rr      = model_sim.extra_params.Fy_rr.data(start_time:end);
Fy_rl      = model_sim.extra_params.Fy_rl.data(start_time:end);
Fy_fr      = model_sim.extra_params.Fy_fr.data(start_time:end);
Fy_fl      = model_sim.extra_params.Fy_fl.data(start_time:end);
Mz_rr      = model_sim.extra_params.Mz_rr.data(start_time:end);
Mz_rl      = model_sim.extra_params.Mz_rl.data(start_time:end);
Mz_fr      = model_sim.extra_params.Mz_fr.data(start_time:end);
Mz_fl      = model_sim.extra_params.Mz_fl.data(start_time:end);
gamma_rr   = model_sim.extra_params.gamma_rr.data(start_time:end);
gamma_rl   = model_sim.extra_params.gamma_rl.data(start_time:end);
gamma_fr   = model_sim.extra_params.gamma_fr.data(start_time:end);
gamma_fl   = model_sim.extra_params.gamma_fl.data(start_time:end);
delta_fr   = model_sim.extra_params.delta_fr.data(start_time:end);
delta_fl   = model_sim.extra_params.delta_fl.data(start_time:end);
mu_f       = model_sim.extra_params.mu_f.data(start_time:end);
mu_r       = model_sim.extra_params.mu_r.data(start_time:end);
Fy__f       = model_sim.extra_params.Fy__f.data(start_time:end);
Fy__r       = model_sim.extra_params.Fy__r.data(start_time:end);
time_sim = time_sim(start_time:end);

% Chassis side slip angle beta [rad]
beta = atan(v./u);

% -----------------
% Accelerations
% -----------------
% Derivatives of u, v [m/s^2]
dot_u = diff(u)/Ts;
dot_v = diff(v)/Ts;
% Total longitudinal and lateral accelerations
Ax = dot_u(1:end) - Omega(2:end).*v(2:end);
Ay = dot_v(1:end) + Omega(2:end).*u(2:end);
Ay_norm = Ay/g;
% Ax low-pass filtered signal (zero-phase digital low-pass filtering)
Wn_filter = 0.01;
[b_butt,a_butt] = butter(4,Wn_filter,'low');
Ax_filt = filtfilt(b_butt,a_butt,Ax);  
dot_u_filt = filtfilt(b_butt,a_butt,dot_u);  
% Steady state lateral acceleration
Ay_ss = Omega.*u;
% Longitudinal jerk [m/s^3]
jerk_x = diff(dot_u)/Ts;

% -----------------
% Other parameters
% -----------------
% Total CoM speed [m/s]
vG = sqrt(u.^2 + v.^2);
% Steady state and transient curvature [m]
rho_ss   = Omega./vG;
rho_tran = ((dot_v.*u(1:end-1) - dot_u.*v(1:end-1)) ./ ((vG(1:end-1)).^3)) + rho_ss(1:end-1);
rho = Omega./u;
% Desired sinusoidal steering angle for the equivalent single track front wheel
desired_steer_atWheel = delta/tau_D;
n_sim = length(model_sim);
FAxc      = CAx*u.^2; % long friction
FAz_f     = CAzf*u.^2;
FAz_r     = CAzr*u.^2;

% LATERAL LOAD TRANSFER
% -----------------
% Longitudinal
DFx_f = (Fx_fr - Fx_fl) / 2;
DFx_r = (Fx_rr - Fx_rl) / 2;

% -----------------
% Lateral
DFy_f = (Fy_fr + Fy_fl) / 2;
DFy_r = (Fy_rr + Fy_rl) / 2;

% -----------------
% Vertical
DFz_f = (Fz_fr - Fz_fl) / 2;
DFz_r = (Fz_rr - Fz_rl) / 2;

% -----------------
% Vertical load
Fzf = m * g * Lr / L - m * Ax * hg / L + FAz_f(1:end-1); % front vertical load
Fzr = m * g * Lf / L + m * Ax * hg / L + FAz_r(1:end-1); % rear vertical load

% NORMALIZED AXLE CHARACTERISTICS
ay_0 = rho_ss .* u; % Fixed acceleration (to be fixed)
ay_0_norm = ay_0 / g;
% -----------------
% alpha rear (rear axle side slio angle)
alpha_r = - (v - Omega * Lr) ./ u;

% normalized rear lateral force
Fz_r0 = m * g * Lf / L;

% -----------------
% alpha front
alpha_f = delta - (v + Omega * Lf) ./ u; % [rad]

% normalized front lateral force
Fz_f0 = m * g * Lr / L;


% Handling diagram
minmax_norm_axle_char = min(max(mu_r), max(mu_f));
maxmin_norm_axle_char = max(min(mu_r), min(mu_f))+2e-3;
Ay_hand = maxmin_norm_axle_char:0.005:minmax_norm_axle_char; % nomralized acc
pos_mu_r_remap = zeros(length(Ay_hand), 1);
pos_mu_f_remap = zeros(length(Ay_hand), 1);
pos_Ay_remap = zeros(length(Ay_hand), 1);
for i=1:size(Ay_hand, 2)
    [~, pos_mu_r_remap(i)] = min(abs(mu_r - Ay_hand(i)));
    [~, pos_mu_f_remap(i)] = min(abs(mu_f - Ay_hand(i)));
    [~, pos_Ay_remap(i)] = min(abs(Ay_norm - Ay_hand(i)));
end
Delta_alpha = alpha_r(pos_mu_r_remap) - alpha_f(pos_mu_f_remap); % [rad]
handling = -Delta_alpha; % [rad]
% handling2 = delta - rho*L;

% Delta_alpha = alpha_r - alpha_f; % [rad]
% rho = (delta_D/tau_D*pi/180 + Delta_alpha)/L; % curvature [1/m]
% [~, idx] = sort(Ay, 'ascend'); % sort the acceleration
% Ay_norm = Ay(idx)/g; % sort and normalize the acceleration
% Delta_alpha = Delta_alpha(idx); % sort the slidign according to the acceleration
% rho = rho(idx); % sort the curvature according to the acceleration
% rho_ss_sort = rho_ss(idx);

% handling = delta_D(1:end-1)/tau_D*pi/180 - rho_ss_sort*L; % [rad] handling metric

% ---------------------------------
%% UNDERSTEERING GRADIENT
% ---------------------------------
% Cornering stiffnesses
Ky_r = diff(Fy__r(pos_mu_r_remap)) ./ diff(alpha_r(pos_mu_r_remap));
Ky_f = diff(Fy__f(pos_mu_f_remap)) ./ diff(alpha_f(pos_mu_f_remap));
C_r = diff(mu_r) ./ diff(alpha_r);
C_f = diff(mu_f) ./ diff(alpha_f);

if sim_options.test_type == 1 % constant velocity
    K_US_theo = - m*g/L*(Lf./Ky_r - Lr./Ky_f);
    K_US_theo2 = - diff(Delta_alpha) ./ diff(Ay_hand');
    K_US_theo3 = diff(delta(1:end-1)*tau_D)./diff(Ay_norm) - L./u(1:end-2).^2;
elseif sim_options.test_type == 2 % constant curvature
    K_US_theo = - m*g/L^2*(Lf./Ky_r - Lr./Ky_f);
    K_US_theo2 = - 1/L*diff(Delta_alpha) ./ diff(Ay_hand');
end
% fitting the linear part
if sim_options.test_type == 1 % constat speed
    Ay_lin_lim = 0.6; % norm. acc. linear limit
elseif sim_options.test_type == 2
    Ay_lin_lim = 0.6;
end
[~, pos_Ay_lin_lim] = min(abs(Ay - g*Ay_lin_lim));
u_lin_lim = u(pos_Ay_lin_lim); % velocity at which we lost the linearity
if max(Ay_hand) < Ay_lin_lim
  ay_max_lin = max(Ay_hand); % [-] maximum norm acceleration for whom linearity holds
else
  ay_max_lin = Ay_lin_lim; % [-] maximum norm acceleration for whom linearity holds 
end
ay_fit = Ay_hand(Ay_hand < ay_max_lin);
ay_fit_lin = [maxmin_norm_axle_char:1e-6:ay_max_lin];
handling_fit = handling(Ay_hand < ay_max_lin);
p_lin = polyfit(ay_fit,handling_fit,1); % fitting the linear part
K_US = p_lin(1); % understeering gradient
handling_fit_lin = polyval(p_lin, ay_fit_lin); % fitted handling diagram
% fitting the nonlinear part 
ay_fit = Ay_hand(Ay_hand > ay_max_lin);
handling_fit = handling(Ay_hand > ay_max_lin) - K_US * ay_fit';
poly_deg = 3;
A = ones(length(ay_fit), poly_deg);
for i = 1:poly_deg - 1 
A(:,i) = ay_fit.^(poly_deg - i + 1);
end
p_tmp = inv(A' * A) * A'*handling_fit; % coefficients of the regression

if ay_max_lin < max(Ay_hand)
  ay_fit_nonlin = [ay_max_lin:0.001:Ay_hand(end)];
  p_nl = p_tmp;
  p_nl(end + 1) = p_nl(end);
  p_nl(end - 1) = K_US; % add the linear term from privious fitting
  handling_fit_nonlin = polyval(p_nl, ay_fit_nonlin); % fitted handling diagram
else
  handling_fit_nonlin = 0;
  ay_fit_nonlin = 0;
end

%----------------------------------
%% CRITICAL SPEED
%----------------------------------
u_cr = sqrt(Ay./rho(1:end-1)); % critical velocity
u_cr_hand = sqrt(Ay_hand*g./rho(pos_Ay_remap)'); % critical velocity at ay resampled
u_slide = sqrt(1/abs(K_US/g));

% ---------------------------------
%% YAW RATE & BETA GAIN
% ---------------------------------
yaw_rate_gain = Omega ./ (delta); % [1/s]
beta_gain = beta ./ (delta); % [-]
yaw_rate_gain_theo = u / L ./ (1 + u.^2 * K_US/g ); % [1/s] theoretical yaw rate gain
u_tmp = u(pos_Ay_remap);
u_tmp = u_tmp(1:end-1);
beta_gain_theo = Lr ./ L - m / L^3 * (Lf^2 ./ Ky_r + Lr^2 ./ Ky_f) .* u_tmp.^2 ./ (1 + K_US/g * u_tmp.^2); % [-] theoretical beta gain

if enable_plot
% ---------------------------------
%% Plot load transfer
% ---------------------------------
% --- DeltaFxf DeltaFxr -- %
  fig_load_transf = figure('Name','Load transfer','NumberTitle','off', 'Color', 'w'); clf
  ax(1) = subplot(221);
  plot(Ay_norm, DFx_f(1:end-1),'LineWidth',2)
  hold on
  plot(Ay_norm, DFx_r(1:end-1), '--', 'LineWidth',2)
  legend('$\Delta F_{xf}$','$\Delta F_{xr}$','location','northwest')
  title('$\Delta F_{xf}$ and $\Delta F_{xr}$ [N]')
  xlabel('$a_y/g [-]$')
  grid on; box on;
  % --- DeltaFyf DeltaFyr -- %
  ax(2) = subplot(222);
  plot(Ay_norm, DFy_f(1:end-1),'LineWidth',2)
  hold on
  plot(Ay_norm, DFy_r(1:end-1), '--', 'LineWidth',2)
  legend('$\Delta F_{yf}$','$\Delta F_{yr}$','location','northwest')
  title('$\Delta F_{yf}$ and $\Delta F_{yr}$ [N]')
  xlabel('$a_y/g [-]$')
  grid on; box on;
  % --- DeltaFzf DeltaFzr -- %
  ax(3) = subplot(223);
  plot(Ay_norm, DFz_f(1:end-1),'LineWidth',2)
  hold on
  plot(Ay_norm, DFz_r(1:end-1), '--', 'LineWidth',2)
  legend('$\Delta F_{zf}$','$\Delta F_{zr}$','location','northwest')
  title('$\Delta F_{zf}$ and $\Delta F_{zr}$ [N]')
  xlabel('$a_y/g [-]$')
  grid on; box on;
  sgtitle('Lateral load transfer', 'FontSize', 25)
  if enable_export == 1
    export_figure(fig_load_transf, '\fig_load_transf.eps', 'images\');
  end
  clear ax

  % % figure('Name', 'Vertical load', 'NumberTitle', 'off'), clf
  % % ax(1) = subplot(1,2,1);
  % % hold on
  % % grid on; box on;
  % % plot(time_sim(1:end-1), Fzf, 'LineWidth', 2) % computed with formulas
  % % plot(time_sim, Fz_fr + Fz_fl, '--', 'LineWidth', 2) % from simulation
  % % title('$F_{zf}$ [N]')
  % % legend('Theo', 'Sim.', 'location', 'northwest')
  % % ax(2) = subplot(1,2,2);
  % % hold on
  % % grid on; box on;
  % % plot(time_sim(1:end-1), Fzr, 'LineWidth', 2) % computed with formulas
  % % plot(time_sim, Fz_rr + Fz_rl, '--', 'LineWidth', 2) % from simulation
  % % title('$F_{zr}$ [N]')
  % % legend('Theo', 'Sim.', 'location', 'northeast')
  % % sgtitle('Vertical load', 'FontSize', 20)

  % % clear ax

  % % ---------------------------------
  % %% Plot axle characteristics
  % % ---------------------------------
  % figure('Name','Long axle char','NumberTitle','off'), clf
  % subplot(1,2,1)
  % hold on
  % grid on; box on;
  % plot(alpha_r, Fx_rr, 'LineWidth',2)
  % plot(alpha_r, Fx_rl, 'LineWidth',2)
  % title('$F_{xr}$ [N]')
  % xlabel('$\alpha_r$')
  % ylabel('$F_{xrr}, F_{xrl}$')
  % legend('$F_{xrr}$','$F_{xrl}$','location','southeast')
  % xlim([0.001 0.06])
  % subplot(1,2,2)
  % hold on
  % grid on; box on;
  % plot(alpha_f, Fx_fr, 'LineWidth',2)
  % plot(alpha_f, Fx_fl, 'LineWidth',2)
  % title('$F_{xf}$ [N]')
  % xlabel('$\alpha_f$')
  % ylabel('$F_{xfr}, F_{xfl}$')
  % legend('$F_{xfr}$','$F_{xfl}$','location','southeast')
  % xlim([0.001 0.06])

  % figure('Name','Lat axle char','NumberTitle','off'), clf
  % subplot(1,2,1)
  % hold on
  % grid on; box on;
  % plot(alpha_r, Fy_rr/Fz_r0, 'LineWidth',2)
  % plot(alpha_r, Fy_rl/Fz_r0, 'LineWidth',2)
  % title('$F_{yr}$ [N]')
  % xlabel('$\alpha_r$')
  % ylabel('$F_{yrr}, F_{yrl}$')
  % legend('$F_{yrr}$','$F_{yrl}$','location','southeast')
  % % xlim([0.001 0.06])
  % subplot(1,2,2)
  % hold on
  % grid on; box on;
  % plot(alpha_f, Fy_fr/Fz_f0, 'LineWidth',2)
  % plot(alpha_f, Fy_fl/Fz_f0, 'LineWidth',2)
  % title('$F_{yf}$ [N]')
  % xlabel('$\alpha_f$')
  % ylabel('$F_{yfr}, F_{yfl}$')
  % legend('$F_{yfr}$','$F_{yfl}$','location','southeast')
  % % xlim([0.001 0.06])

  % ---------------------------------
  %% Plot normalized axle characteristics
  % ---------------------------------
  % figure('Name','Axle char. ','NumberTitle','off'), clf
  % hold on
  % grid on; box on;
  % plot(alpha_r, Fy__r, 'LineWidth',2)
  % plot(alpha_f, Fy__f, 'LineWidth',2)
  % title('$\Fy_r, \Fy_f $')
  % xlabel('$\alpha_r, \alpha_f [deg]$')
  % ylabel('[N]')
  % legend('$\Fy_r$', '$\Fy_f $','location','best')

  % --- mu_r -- %
  % idx = time_sim > 2;
  fig_norm_axle_char = figure('Name','Norm axle char. 2','NumberTitle','off', 'Color', 'w'); clf
  hold on; grid on; box on;
  plot(alpha_r*180/pi, mu_r, 'LineWidth',2)
  plot(alpha_f*180/pi, mu_f, 'LineWidth',2)
  title('$\mu_r, \mu_f $')
  xlabel('$\alpha_r, \alpha_f [deg]$')
  ylabel('$\mu_r, \mu_f$ [-]')
  legend('$\mu_r$','$\mu_f$','location','best')
  if enable_export == 1
    export_figure(fig_norm_axle_char, '\fig_norm_axle_char.eps', 'images\');
  end


  % ---------------------------------
  %% Plot handling digram
  % ---------------------------------
  fig_hand = figure('Name','Handling diagram','NumberTitle','off', 'Color', 'w'); clf
  hold on
  grid on; box on;
  plot(Ay_hand, handling, 'LineWidth',2, 'DisplayName','Experimental')
%   plot(Ay_norm, handling2(1:end-1), 'LineWidth',2, 'DisplayName','Handling 2')
  %plot(Ay_hand, -Delta_alpha, 'LineWidth',2)
  plot(ay_fit_lin, handling_fit_lin, '--', 'LineWidth',2.5, 'Displayname','Fit in linear range')
  plot(ay_fit_nonlin, handling_fit_nonlin, '--', 'LineWidth',2.5, 'Displayname','Fit in non-linear range')
  title('Handling diagram')
  xlabel('$a_{y}/g$ [m/s$^2$]')
  ylabel('$\delta_{D}\tau_{H} - \rho L \ [rad]$')
  legend('location', 'northeast')
  if enable_export == 1
    export_figure(fig_hand, '\fig_hand.eps', 'images\');
  end

  % ---------------------------------
  %% Plot understeering gradient
  % ---------------------------------
  fig_KUS = figure('Name','Understeering grad','NumberTitle','off', 'Color','w'); clf
  hold on
  grid on; box on;
  if sim_options.test_type == 1
    plot(Ay_hand(1:end-1), K_US_theo, 'LineWidth',2, 'DisplayName','$-\frac{mg}{L}(\frac{L_f}{K_{yr}} - \frac{Lr}{K_{yf}})$')
    plot(Ay_hand(1:end-1), K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$- \frac{d\Delta \alpha}{d ay/g}$')
    plot(ay_fit_lin, K_US*ones(length(ay_fit_lin), 1), '--r', 'LineWidth',2, 'DisplayName', 'Linear $K_{US}$');
  elseif sim_options.test_type == 2
    plot(Ay_hand(1:end-1), K_US_theo, 'LineWidth',2, 'DisplayName','$-\frac{mg}{L^2}(\frac{L_f}{K_{yr}} - \frac{Lr}{K_{yf}})$')
    plot(Ay_hand(1:end-1), K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$- \frac{1}{L}\frac{d\Delta \alpha}{d ay/g}$')
    plot(ay_fit_lin, K_US*ones(length(ay_fit_lin), 1)/L, '--r', 'LineWidth',2, 'DisplayName', 'Linear $K_{US}$');
  end
  %plot(Ay_norm(1:end-1), K_US_theo3, 'LineWidth',2, 'DisplayName','Formula 2')
  title('Normalized understeering gradient')
  legend('location', 'northeast')
  xlim("padded")
  ylabel('g$K_{US}$')
  xlabel('ay/g [-]')
  if enable_export == 1
    export_figure(fig_KUS, '\fig_KUS.eps', 'images\');
  end

  % ---------------------------------
  %% PLOT BETA AND YAW RATE GAINS
  % ---------------------------------
  fig_yaw_gain = figure('Name','Yaw rate gain','NumberTitle','off', 'Color', 'w'); clf
  hold on 
  grid on; box on;
  plot(u*3.6, yaw_rate_gain, 'LineWidth',2, 'DisplayName','$\frac{\Omega}{\delta}$')
  plot(u*3.6, yaw_rate_gain_theo, 'LineWidth',2, 'DisplayName','$\frac{u}{L(1+K_{US}u^2)}$')
  xline(u_lin_lim*3.6, '--r', 'LineWidth',2, 'DisplayName','$K_{US}$ lin limit')
  title('Yaw rate gain')
  xlabel('u [km/h]')
  ylabel('$\frac{\Omega}{\delta}$')
  legend('location', 'northwest')
  if enable_export == 1
    export_figure(fig_yaw_gain, '\fig_yaw_gain.eps', 'images\');
  end

  fig_beta_gain = figure('Name','Beta gain','NumberTitle','off'); clf
  hold on 
  grid on; box on;
  plot(u*3.6, beta_gain, 'LineWidth',2, 'DisplayName','$\frac{\beta}{\delta}$')
  plot(u_tmp*3.6, beta_gain_theo, 'LineWidth',2, 'DisplayName','Theo.')
  xline(u_lin_lim*3.6, 'r--', 'LineWidth',2, 'DisplayName','$K_{US}$ linear limit')
  title('Body slip angle gain')
  xlabel('u [km/h]')
  ylabel('$\frac{\beta}{\delta}$ ')
  legend('location', 'southwest')
  if enable_export == 1
    export_figure(fig_beta_gain, '\fig_beta_gain.eps', 'images\');
  end

  % % % ---------------------------------
  % % %% HANDLING AS FUNCTION OF THE TIME
  % % % ---------------------------------
  % figure()
  % hold on
  % plot(time_sim, mu_r, 'DisplayName','$\mu_r$')
  % plot(time_sim, mu_f, 'DisplayName','$\mu_f$')
  % legend('Location','southeast')
  % xlabel('Time [s]')
  % ylabel('$\mu_r$, $\mu_f$')
  % grid on; box on;
end
%% SAVE DATA ON FILE
data.Ay_hand = Ay_hand;
data.handling = handling;
data.K_US_theo2 = K_US_theo2;
save(strcat("saved_test\", name_output), 'data');

