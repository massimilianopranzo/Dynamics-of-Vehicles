function dataAnalysisAxleChar(model_sim,vehicle_data,Ts)
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


    %% Print on the terminal
    fprintf("mu_r = %.4f\n", mu_r(end));
    fprintf("mu_f = %.4f\n", mu_f(end));
    fprintf("alfa_f = %.4f\n", alpha_f(end));
    fprintf("alfa_r = %.4f\n", alpha_r(end));
    fprintf("ay = %.4f\n", Ay(end));
pause 
    %% Print on file
    fileID = fopen('results_axle_char.txt','a');
    % pedal, steering angle, mu_r, mu_f, alfa_f, alfa_r, ay
fprintf(fileID,'%.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n', ped_0(end), delta_D(end), mu_r(end), mu_f(end), alpha_f(end), alpha_r(end), Ay(end));
fclose(fileID);
    
end
    