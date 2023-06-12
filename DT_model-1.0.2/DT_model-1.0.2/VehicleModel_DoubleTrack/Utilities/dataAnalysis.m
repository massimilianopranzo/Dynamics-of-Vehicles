function dataAnalysis(model_sim,vehicle_data,Ts)
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

    % LATERAL LOAD TRANSFER
    % -----------------
    % Longitudinal load transfer
    DFx_f = (Fx_fr - Fx_fl) / 2;
    DFx_r = (Fx_rr - Fx_rl) / 2;

    % -----------------
    % Lateral load transfer
    DFy_f = (Fy_fr + Fy_fl) / 2;
    DFy_r = (Fy_rr + Fy_rl) / 2;

    % -----------------
    % Vertical load transfer
    DFz_f = (Fz_fr - Fz_fl) / 2;
    DFz_r = (Fz_rr - Fz_rl) / 2;

    % -----------------
    % Vertical load
    Fzf = m * g * Lr / L - m * Ax * hg / L + FAz_f(1:end-1); % front vertical load
    Fzr = m * g * Lf / L + m * Ax * hg / L + FAz_r(1:end-1); % rear vertical load

    % -----------------
    % Lateral load


    % NORMALIZED AXLE CHARACTERISTICS
    ay_0 = rho_ss .* u; % Fixed acceleration (to be fixed)
    ay_0_norm = ay_0 / g;
    % -----------------
    % alpha rear (rear axle side slio angle)
    alpha_r = - (v - Omega * Lr) ./ u;

    % normalized rear lateral force
    Fz_r0 = m * g * Lf / L;
    Fy_r = Fy_rr + Fy_rl;
    mu_r = Fy_r / Fz_r0;

    % -----------------
    % alpha front
    alpha_f = delta_D/tau_D*pi/180 - (v + Omega * Lf) ./ u;

    % normalized front lateral force
    Fz_f0 = m * g * Lr / L;
    Fy_f = sin(delta_fl).*Fx_fl + cos(delta_fl).*Fy_fl + sin(delta_fr).*Fx_fr + cos(delta_fr).*Fy_fr;
    mu_f = Fy_f / Fz_f0;


    % Handling diagram
    Delta_alpha = alpha_r - alpha_f; % [rad]
    rho = (delta_D/tau_D*pi/180 + Delta_alpha)/L; % curvature [1/m]
    [~, idx] = sort(Ay, 'ascend'); % sort the acceleration
    Ay_norm = Ay(idx)/g; % sort and normalize the acceleration
    Delta_alpha = Delta_alpha(idx); % sort the slidign according to the acceleration
    rho = rho(idx); % sort the curvature according to the acceleration
    rho_ss_sort = rho_ss(idx);

    handling = delta_D(1:end-1)/tau_D*pi/180 - rho*L; % [rad] handling metric
    
    % ---------------------------------
    %% UNDERSTEERING GRADIENT
    % ---------------------------------
    % Cornering stiffnesses
    C_alpha_r = diff(Fy_r) ./ diff(alpha_r);
    C_alpha_f = diff(Fy_f) ./ diff(alpha_f);
    C_alpha_r = C_alpha_r(idx);
    C_alpha_f = C_alpha_f(idx);

    K_US_theo = - m / (L / tau_D) * (Lf ./ C_alpha_r - Lr ./ C_alpha_f);
    K_US_theo2 = - diff(Delta_alpha) ./ diff(Ay);

    % fitting the linear part
    ay_max_lin = 1e-4; % [-] maximum norm acceleration for whom linearity holds 
    ay_fit = Ay_norm(Ay_norm < ay_max_lin);
    ay_fit_lin = [0:1e-6:ay_max_lin];
    handling_fit = handling(Ay_norm < ay_max_lin);
    p_lin = polyfit(ay_fit,handling_fit,1); % fitting the linear part
    K_US = p_lin(1); % understeering gradient
    handling_fit_lin = polyval(p_lin, ay_fit_lin); % fitted handling diagram
    % fitting the nonlinear part 
    ay_fit = Ay_norm(Ay_norm > ay_max_lin);
    handling_fit = handling(Ay_norm > ay_max_lin) - K_US * ay_fit;
    poly_deg = 3;
    A = ones(length(ay_fit), poly_deg);
    for i = 1:poly_deg - 1 
        A(:,i) = ay_fit.^(poly_deg - i + 1);
    end
    p_tmp = inv(A' * A) * A'*handling_fit; % coefficients of the regression
    
    ay_fit_nonlin = [ay_max_lin:0.001:Ay_norm(end)];
    p_nl = p_tmp;
    p_nl(end + 1) = p_nl(end);
    p_nl(end - 1) = K_US; % add the linear term from privious fitting
    handling_fit_nonlin = polyval(p_nl, ay_fit_nonlin); % fitted handling diagram

    % ---------------------------------
    %% YAW RATE & BETA GAIN
    % ---------------------------------
    yaw_rate_gain = Omega ./ (delta_D * pi / 180); % [1/s]
    beta_gain = beta ./ (delta_D*pi/180); % [-]
    yaw_rate_gain_theo = u / L / tau_D ./ (1 + u.^2 * K_US); % [1/s] theoretical yaw rate gain
    beta_gain_theo = Lr ./ (tau_D * L) - m / L^3 * (Lf^2 ./ C_alpha_r + Lr^2 ./ C_alpha_f) ./ tau_D .* (u(1:end-1).^2 ./ (1 + K_US * u(1:end-1).^2)); % [-] theoretical beta gain



    
    % ---------------------------------
    %% PLOTS
    % ---------------------------------

    % ---------------------------------
    %% Plot vehicle inputs
    % ---------------------------------
    figure('Name','Inputs','NumberTitle','off'), clf   
    % --- pedal --- %
    ax(1) = subplot(211);
    hold on
    plot(time_sim,ped_0,'LineWidth',2)
    grid on
    title('pedal $p_0$ [-]')
    xlim([0 time_sim(end)])
    % --- delta_0 --- %
    ax(2) = subplot(212);
    plot(time_sim,delta_D,'LineWidth',2)
    grid on
    title('steering angle $\delta_D$ [deg]')
    xlim([0 time_sim(end)])
    
    % ---------------------------------
    %% Plot vehicle motion
    % ---------------------------------
    figure('Name','Motion','NumberTitle','off'), clf   
    % --- u --- %
    ax(1) = subplot(221);
    plot(time_sim,u*3.6,'LineWidth',2)
    grid on
    title('$u$ [km/h]')
    xlim([0 time_sim(end)])
    % --- v --- %
    ax(2) = subplot(222);
    plot(time_sim,v,'LineWidth',2)
    grid on
    title('$v$ [m/s]')
    xlim([0 time_sim(end)])
    % --- Omega --- %
    ax(3) = subplot(223);
    plot(time_sim,Omega,'LineWidth',2)
    grid on
    title('$\Omega$ [rad/s]')
    xlim([0 time_sim(end)])

    % % ---------------------------------
    % %% Plot steering angles
    % % ---------------------------------
    % figure('Name','Steer','NumberTitle','off'), clf   
    % % --- delta_0 --- %
    % ax(1) = subplot(221);
    % plot(time_sim,delta_D,'LineWidth',2)
    % grid on
    % title('$\delta_0$ [deg]')
    % xlim([0 time_sim(end)])
    % % --- delta_fr --- %
    % ax(2) = subplot(222);
    % plot(time_sim,delta_fr,'LineWidth',2)
    % grid on
    % title('$\delta_{fr}$ [deg]')
    % xlim([0 time_sim(end)])
    % % --- delta_fl --- %
    % ax(3) = subplot(223);
    % hold on
    % plot(time_sim,delta_fl,'LineWidth',2)
    % grid on
    % title('$\delta_{fl}$ [deg]')
    % xlim([0 time_sim(end)])
    % % --- comparison --- %
    % ax(4) = subplot(224);
    % hold on
    % plot(time_sim,delta_D/tau_D,'LineWidth',2)
    % plot(time_sim,delta_fr,'LineWidth',2)
    % plot(time_sim,delta_fl,'LineWidth',2)
    % grid on
    % legend('$\delta_D/\tau_D$','$\delta_{fr}$','$\delta_{fl}$','location','best')
    % xlim([0 time_sim(end)])

    % -------------------------------
    %% Plot lateral tire slips and lateral forces
    % -------------------------------
    figure('Name','Lat slips & forces','NumberTitle','off'), clf
    % --- alpha_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,alpha_rr,'LineWidth',2)
    grid on
    title('$\alpha_{rr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,alpha_rl,'LineWidth',2)
    grid on
    title('$\alpha_{rl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,alpha_fr,'LineWidth',2)
    grid on
    title('$\alpha_{fr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,alpha_fl,'LineWidth',2)
    grid on
    title('$\alpha_{fl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- Fy_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Fy_rr,'LineWidth',2)
    grid on
    title('$Fy_{rr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fy_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Fy_rl,'LineWidth',2)
    grid on
    title('$Fy_{rl}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Fy_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Fy_fr,'LineWidth',2)
    grid on
    title('$Fy_{fr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fy_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Fy_fl,'LineWidth',2)
    grid on
    title('$Fy_{fl}$ [N]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    
    % % ---------------------------------
    % %% Plot longitudinal tire slips and longitudinal forces
    % % ---------------------------------
    % figure('Name','Long slips & forces','NumberTitle','off'), clf
    % % --- kappa_rr --- %
    % ax(1) = subplot(331);
    % plot(time_sim,kappa_rr,'LineWidth',2)
    % grid on
    % title('$\kappa_{rr}$ [-]')
    % xlim([0 time_sim(end)])
    % % --- kappa_rl --- %
    % ax(2) = subplot(332);
    % plot(time_sim,kappa_rl,'LineWidth',2)
    % grid on
    % title('$\kappa_{rl}$ [-]')
    % xlim([0 time_sim(end)])
    % % --- kappa_fr --- %
    % ax(3) = subplot(333);
    % plot(time_sim,kappa_fr,'LineWidth',2)
    % grid on
    % title('$\kappa_{fr}$ [-]')
    % xlim([0 time_sim(end)])
    % % --- kappa_fl --- %
    % ax(4) = subplot(334);
    % plot(time_sim,kappa_fl,'LineWidth',2)
    % grid on
    % title('$\kappa_{fl}$ [-]')
    % xlim([0 time_sim(end)])
    % % --- Fx_rr --- %
    % ax(5) = subplot(335);
    % plot(time_sim,Fx_rr,'LineWidth',2)
    % grid on
    % title('$Fx_{rr}$ [N]')
    % xlim([0 time_sim(end)])
    % % --- Fx_rl --- %
    % ax(6) = subplot(336);
    % plot(time_sim,Fx_rl,'LineWidth',2)
    % grid on
    % title('$Fx_{rl}$ [N]')
    % xlim([0 time_sim(end)])
    % % --- Fx_fr --- %
    % ax(7) = subplot(337);
    % plot(time_sim,Fx_fr,'LineWidth',2)
    % grid on
    % title('$Fx_{fr}$ [N]')
    % xlim([0 time_sim(end)])
    % % --- Fx_fl --- %
    % ax(8) = subplot(338);
    % plot(time_sim,Fx_fl,'LineWidth',2)
    % grid on
    % title('$Fx_{fl}$ [N]')
    % xlim([0 time_sim(end)])
    
    % % linkaxes(ax,'x')
    % clear ax

    % % ---------------------------------
    % %% Plot wheel torques and wheel rates
    % % ---------------------------------
    % figure('Name','Wheel rates & torques','NumberTitle','off'), clf
    % % --- omega_rr --- %
    % ax(1) = subplot(331);
    % plot(time_sim,omega_rr,'LineWidth',2)
    % grid on
    % title('$\omega_{rr}$ [rad/s]')
    % xlim([0 time_sim(end)])
    % % --- omega_rl --- %
    % ax(2) = subplot(332);
    % plot(time_sim,omega_rl,'LineWidth',2)
    % grid on
    % title('$\omega_{rl}$ [rad/s]')
    % xlim([0 time_sim(end)])
    % % --- omega_fr --- %
    % ax(3) = subplot(333);
    % plot(time_sim,omega_fr,'LineWidth',2)
    % grid on
    % title('$\omega_{fr}$ [rad/s]')
    % xlim([0 time_sim(end)])
    % % --- omega_fl --- %
    % ax(4) = subplot(334);
    % plot(time_sim,omega_fl,'LineWidth',2)
    % grid on
    % title('$\omega_{fl}$ [rad/s]')
    % xlim([0 time_sim(end)])
    % % --- Tw_rr --- %
    % ax(5) = subplot(335);
    % plot(time_sim,Tw_rr,'LineWidth',2)
    % grid on
    % title('$Tw_{rr}$ [Nm]')
    % xlim([0 time_sim(end)])
    % % --- Tw_rl --- %
    % ax(6) = subplot(336);
    % plot(time_sim,Tw_rl,'LineWidth',2)
    % grid on
    % title('$Tw_{rl}$ [Nm]')
    % xlim([0 time_sim(end)])
    % % --- Tw_fr --- %
    % ax(7) = subplot(337);
    % plot(time_sim,Tw_fr,'LineWidth',2)
    % grid on
    % title('$Tw_{fr}$ [Nm]')
    % xlim([0 time_sim(end)])
    % % --- Tw_fl --- %
    % ax(8) = subplot(338);
    % plot(time_sim,Tw_fl,'LineWidth',2)
    % grid on
    % title('$Tw_{fl}$ [Nm]')
    % xlim([0 time_sim(end)])

    % % linkaxes(ax,'x')
    % clear ax

    % % ---------------------------------
    % %% Plot vertical tire loads and self-aligning torques
    % % ---------------------------------
    % figure('Name','Vert loads & ali torques','NumberTitle','off'), clf
    % % --- Fz_rr --- %
    % ax(1) = subplot(331);
    % plot(time_sim,Fz_rr,'LineWidth',2)
    % grid on
    % title('$Fz_{rr}$ [N]')
    % xlim([0 time_sim(end)])
    % % --- Fz_rl --- %
    % ax(2) = subplot(332);
    % plot(time_sim,Fz_rl,'LineWidth',2)
    % grid on
    % title('$Fz_{rl}$ [N]')
    % xlim([0 time_sim(end)])
    % % --- Fz_fr --- %
    % ax(3) = subplot(333);
    % plot(time_sim,Fz_fr,'LineWidth',2)
    % grid on
    % title('$Fz_{fr}$ [N]')
    % xlim([0 time_sim(end)])
    % % --- Fz_fl --- %
    % ax(4) = subplot(334);
    % plot(time_sim,Fz_fl,'LineWidth',2)
    % grid on
    % title('$Fz_{fl}$ [N]')
    % xlim([0 time_sim(end)])
    % % --- Mz_rr --- %
    % ax(5) = subplot(335);
    % plot(time_sim,Mz_rr,'LineWidth',2)
    % grid on
    % title('$Mz_{rr}$ [Nm]')
    % xlim([0 time_sim(end)])
    % % --- Mz_rl --- %
    % ax(6) = subplot(336);
    % plot(time_sim,Mz_rl,'LineWidth',2)
    % grid on
    % title('$Mz_{rl}$ [Nm]')
    % xlim([0 time_sim(end)])
    % % --- Mz_fr --- %
    % ax(7) = subplot(337);
    % plot(time_sim,Mz_fr,'LineWidth',2)
    % grid on
    % title('$Mz_{fr}$ [Nm]')
    % xlim([0 time_sim(end)])
    % % --- Mz_fl --- %
    % ax(8) = subplot(338);
    % plot(time_sim,Mz_fl,'LineWidth',2)
    % grid on
    % title('$Mz_{fl}$ [Nm]')
    % xlim([0 time_sim(end)])

    % % linkaxes(ax,'x')
    % clear ax

    
    % % ---------------------------------
    % %% Plot wheel camber
    % % ---------------------------------
    % figure('Name','Wheel camber','NumberTitle','off'), clf
    % % --- gamma_rr --- %
    % ax(1) = subplot(221);
    % plot(time_sim,gamma_rr,'LineWidth',2)
    % grid on
    % title('$\gamma_{rr}$ [deg]')
    % xlim([0 time_sim(end)])
    % % --- gamma_rl --- %
    % ax(2) = subplot(222);
    % plot(time_sim,gamma_rl,'LineWidth',2)
    % grid on
    % title('$\gamma_{rl}$ [deg]')
    % xlim([0 time_sim(end)])
    % % --- gamma_fr --- %
    % ax(3) = subplot(223);
    % plot(time_sim,gamma_fr,'LineWidth',2)
    % grid on
    % title('$\gamma_{fr}$ [deg]')
    % xlim([0 time_sim(end)])
    % % --- gamma_fl --- %
    % ax(4) = subplot(224);
    % plot(time_sim,gamma_fl,'LineWidth',2)
    % grid on
    % title('$\gamma_{fl}$ [deg]')
    % xlim([0 time_sim(end)])

    % % linkaxes(ax,'x')
    % clear ax
    
    % ---------------------------------
    %% Plot accelerations, chassis side slip angle and curvature
    % ---------------------------------
    figure('Name','Extra','NumberTitle','off'), clf
    % --- ax --- %
    ax(1) = subplot(221);
    plot(time_sim(2:end),dot_u - Omega(2:end).*v(2:end),'LineWidth',2)
    hold on
    plot(time_sim(2:end),diff(u)/Ts,'--g','LineWidth',2)
    plot(time_sim(2:end),Ax_filt,'-.b','LineWidth',1)
    plot(time_sim(2:end),dot_u_filt,'-.r','LineWidth',1)
    grid on
    title('$a_{x}$ $[m/s^2]$')
    legend('$\dot{u}-\Omega v$','$\dot{u}$','filt $\dot{u}-\Omega v$','filt $\dot{u}$','Location','northeast')
    xlim([0 time_sim(end)])
    % --- ay --- %
    ax(2) = subplot(222);
    plot(time_sim(2:end),dot_v + Omega(2:end).*u(2:end),'LineWidth',2)
    hold on
    plot(time_sim(2:end),Omega(2:end).*u(2:end),'--g','LineWidth',1)
    grid on
    title('$a_{y}$ $[m/s^2]$')
    legend('$\dot{v}+\Omega u$','$\Omega u$','Location','best')
    xlim([0 time_sim(end)])
    % --- beta --- %
    ax(3) = subplot(223);
    plot(time_sim,rad2deg(beta),'LineWidth',2)
    grid on
    title('$\beta$ [deg]')
    xlim([0 time_sim(end)])
    % --- rho --- %
    ax(4) = subplot(224);
    plot(time_sim,rho_ss,'LineWidth',2)
    hold on
    plot(time_sim(1:end-1),rho_tran,'--g','LineWidth',1)
    grid on
    title('$\rho$ [$m^{-1}$]')
    legend('$\rho_{ss}$','$\rho_{transient}$','Location','best')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    % % ---------------------------------
    % %% Plot vehicle pose x,y,psi
    % % ---------------------------------
    % figure('Name','Pose','NumberTitle','off'), clf 
    % % --- x --- %
    % ax(1) = subplot(221);
    % plot(time_sim,x_CoM,'LineWidth',2)
    % grid on
    % title('$x$ [m]')
    % xlim([0 time_sim(end)])
    % % --- y --- %
    % ax(2) = subplot(222);
    % plot(time_sim,y_CoM,'LineWidth',2)
    % grid on
    % title('$y$ [m]')
    % xlim([0 time_sim(end)])
    % % --- psi --- %
    % ax(3) = subplot(223);
    % plot(time_sim,rad2deg(psi),'LineWidth',2)
    % grid on
    % title('$\psi$ [deg]')
    % xlim([0 time_sim(end)])

    % % linkaxes(ax,'x')
    % clear ax

    % % -------------------------------
    % %% Plot G-G diagram from simulation data
    % % -------------------------------
    % figure('Name','G-G plot','NumberTitle','off'), clf
    % axis equal
    % hold on
    % plot3(Ay,Ax_filt,u(1:end-1),'Color',color('purple'),'LineWidth',3)
    % xlabel('$a_y$ [m/s$^2$]')
    % ylabel('$a_x$ [m/s$^2$]')
    % zlabel('$u$ [m/s]')
    % title('G-G diagram from simulation data','FontSize',18)
    % grid on

    % -------------------------------
    %% Plot vehicle path
    % -------------------------------
    N = length(time_sim);
    figure('Name','Real Vehicle Path','NumberTitle','off'), clf
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
    grid on
    hold off
    
    % % ---------------------------------
    % %% Plot load transfer
    % % ---------------------------------
    % % --- DeltaFxf DeltaFxr -- %
    % figure('Name','Load transfer','NumberTitle','off'), clf
    % ax(1) = subplot(221);
    % plot(time_sim, DFx_f,'LineWidth',2)
    % hold on
    % plot(time_sim, DFx_r, '--', 'LineWidth',2)
    % legend('$\Delta F_{xf}$','$\Delta F_{xr}$','location','best')
    % title('$\Delta F_{xf}$ and $\Delta F_{xr}$ [N]')
    % grid on
    % % --- DeltaFyf DeltaFyr -- %
    % ax(2) = subplot(222);
    % plot(time_sim, DFy_f,'LineWidth',2)
    % hold on
    % plot(time_sim, DFy_r, '--', 'LineWidth',2)
    % legend('$\Delta F_{yf}$','$\Delta F_{yr}$','location','best')
    % title('$\Delta F_{yf}$ and $\Delta F_{yr}$ [N]')
    % grid on
    % % --- DeltaFzf DeltaFzr -- %
    % ax(3) = subplot(223);
    % plot(time_sim, DFz_f,'LineWidth',2)
    % hold on
    % plot(time_sim, DFz_r, '--', 'LineWidth',2)
    % legend('$\Delta F_{zf}$','$\Delta F_{zr}$','location','best')
    % title('$\Delta F_{zf}$ and $\Delta F_{zr}$ [N]')
    % grid on
    % sgtitle('Lateral load transfer', 'FontSize', 20)
    
    % clear ax

    % % figure('Name', 'Vertical load', 'NumberTitle', 'off'), clf
    % % ax(1) = subplot(1,2,1);
    % % hold on
    % % grid on
    % % plot(time_sim(1:end-1), Fzf, 'LineWidth', 2) % computed with formulas
    % % plot(time_sim, Fz_fr + Fz_fl, '--', 'LineWidth', 2) % from simulation
    % % title('$F_{zf}$ [N]')
    % % legend('Theo', 'Sim.', 'location', 'northwest')
    % % ax(2) = subplot(1,2,2);
    % % hold on
    % % grid on
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
    % grid on
    % plot(alpha_r, Fx_rr, 'LineWidth',2)
    % plot(alpha_r, Fx_rl, 'LineWidth',2)
    % title('$F_{xr}$ [N]')
    % xlabel('$\alpha_r$')
    % ylabel('$F_{xrr}, F_{xrl}$')
    % legend('$F_{xrr}$','$F_{xrl}$','location','southeast')
    % xlim([0.001 0.06])
    % subplot(1,2,2)
    % hold on
    % grid on
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
    % grid on
    % plot(alpha_r, Fy_rr/Fz_r0, 'LineWidth',2)
    % plot(alpha_r, Fy_rl/Fz_r0, 'LineWidth',2)
    % title('$F_{yr}$ [N]')
    % xlabel('$\alpha_r$')
    % ylabel('$F_{yrr}, F_{yrl}$')
    % legend('$F_{yrr}$','$F_{yrl}$','location','southeast')
    % % xlim([0.001 0.06])
    % subplot(1,2,2)
    % hold on
    % grid on
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
    % --- mu_r -- %
    data_read = readmatrix('results_axle_char.txt');
    ped_0_file = data_read(:,1);
    delta_D_file = data_read(:,2);
    mu_r_file = data_read(:,3);
    mu_f_file = data_read(:,4);
    alpfa_f_file = data_read(:,5);
    alpfa_r_file = data_read(:,6);
    Ay_file = data_read(:,7);
    idx = time_sim > 2;
    figure('Name','Norm axle char. 2','NumberTitle','off'), clf
    hold on
    grid on
    plot(alpha_r(idx)*180/pi, mu_r(idx), 'LineWidth',2)
    plot(alpha_f(idx)*180/pi, mu_f(idx), 'LineWidth',2)
    plot(alpfa_r_file*180/pi, mu_r_file,  'o', 'LineWidth', 2)
    plot(alpfa_f_file*180/pi, mu_f_file,  'o', 'LineWidth', 2)
    title('$\mu_r, \mu_f$')
    xlabel('$\alpha_r, \alpha_f$')
    ylabel('$\mu_r, \mu_f$')
    legend('$\mu_r$','$\mu_f$','$\mu_r ss$', '$\mu_f ss$', 'location','best')
    

    % % ---------------------------------
    % %% Plot handling digram
    % % ---------------------------------
    % figure('Name','Handling diagram','NumberTitle','off'), clf
    % hold on
    % grid on
    % plot(Ay_norm, handling, 'LineWidth',2)
    % plot(ay_fit_lin, handling_fit_lin, '--', 'LineWidth',2)
    % plot(ay_fit_nonlin, handling_fit_nonlin, '--', 'LineWidth',2)
    % title('Handling diagram')
    % ylabel('$\delta_{D}\tau_{H} - \rho L [rad]$')
    % legend('Data', 'Fit in linear range', 'Fit in non-linear range', 'location', 'northeast')
    
    % % ---------------------------------
    % %% Plot understeering gradient
    % % ---------------------------------
    % figure('Name','Understeering grad','NumberTitle','off'), clf
    % hold on
    % grid on
    % plot(Ay_norm, K_US_theo, 'LineWidth',2)
    % plot(Ay_norm(1:end-1), K_US_theo2, 'LineWidth',2)
    % title('Understeering gradient')
    % legend('Formula', 'Diff', 'location', 'northeast')
    % xlim("padded")
    % ylabel('$K_{US}$')


    % % ---------------------------------
    % %% PLOT BETA AND YAW RATE GAINS
    % % ---------------------------------
    % figure('Name','Yaw rate gain','NumberTitle','off'), clf
    % hold on 
    % grid on
    % plot(u, yaw_rate_gain, 'LineWidth',2)
    % plot(u, yaw_rate_gain_theo, 'LineWidth',2)
    % plot(u, beta_gain, 'LineWidth',2)
    % plot(u(1:end-1), beta_gain_theo, 'LineWidth',2)
    % title('Yaw rate and $\beta$ gain')
    % legend('Yaw gain measure', 'Yaw gain theo', '$\beta$ gain measure', '$\beta$ gain theo', 'location', 'northeast')

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
    % grid on
    % title('$F_{rr}$ [deg]')
    % % --- alpha_rl --- %
    % ax(2) = subplot(223);
    % plot(alpha_rl,Fy_rl,'LineWidth',2)
    % grid on
    % title('$F_{rl}$ [deg]')
    % % --- alpha_fr --- %
    % ax(3) = subplot(222);
    % plot(alpha_fr,Fy_fr,'LineWidth',2)
    % grid on
    % title('$F_{fr}$ [deg]')
    % % --- alpha_fl --- %
    % ax(4) = subplot(221);
    % plot(alpha_fl,Fy_fl,'LineWidth',2)
    % grid on
    % title('$F_{fl}$ [deg]')
    % sgtitle('Lateral load $F_{y}$')
    % % linkaxes(ax,'x')
    % clear ax
end
    