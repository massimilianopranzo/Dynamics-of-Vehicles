function analysis_store = dataAnalysisMulti(model_sim, vehicle_data, Ts)
    close all
    % ----------------------------------------------------------------
    %% Post-Processing and Data Analysis
    % ----------------------------------------------------------------

    % ---------------------------------
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
    delta_f0 = vehicle_data.front_wheel.delta_f0;  % [-] front toe angle
    delta_r0 = vehicle_data.rear_wheel.delta_r0;  % [-] rear toe angle
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

    % LATERAL LOAD TRANSFER
    % -----------------
    % Longitudinal load transfer
    DFx_f = (Fx_fr - Fx_fl) / 2;
    DFx_r = (Fx_rr - Fx_rl) / 2;

    % -----------------
    % Lateral load transfer
    DFy_f = (Fy_fr - Fy_fl) / 2;
    DFy_r = (Fy_rr - Fy_rl) / 2;

    % -----------------
    % Vertical load transfer
    DFz_f = (Fz_fr - Fz_fl) / 2;
    DFz_r = (Fz_rr - Fz_rl) / 2;

    % NORMALIZED AXLE CHARACTERISTICS
    ay_0 = rho_ss .* u; % Fixed acceleration (to be fixed)
    ay_0_norm = ay_0 / g;
    % -----------------
    % alpha rear (rear axle side slio angle)
    alpha_r = delta_r0 - (v - Omega * Lr) ./ u;

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
    %% Plot steering angles
    % ---------------------------------
    analysis_store.time_sim 			= time_sim;
    analysis_store.delta_D 				= delta_D;
    analysis_store.delta_fr 			= delta_fr;
    analysis_store.delta_fl 			= delta_fl;
    analysis_store.alpha_rr 			= alpha_rr;
    analysis_store.alpha_rl 			= alpha_rl;
    analysis_store.alpha_fl 			= alpha_fl;
    analysis_store.alpha_fr 			= alpha_fr;
    analysis_store.alpha_r 				= alpha_r;
    analysis_store.alpha_f 				= alpha_f;
    analysis_store.Fy_rr 				= Fy_rr;
    analysis_store.Fy_rl 				= Fy_rl;
    analysis_store.Fy_fl 				= Fy_fl;
    analysis_store.Fy_fr 				= Fy_fr;
    analysis_store.kappa_rr 			= kappa_rr;
    analysis_store.kappa_rl 			= kappa_rl;
    analysis_store.kappa_fl 			= kappa_fl;
    analysis_store.kappa_fr 			= kappa_fr;
    analysis_store.Fx_rr 				= Fx_rr;
    analysis_store.Fx_rl 				= Fx_rl;
    analysis_store.Fx_fl 				= Fx_fl;
    analysis_store.Fx_fr 				= Fx_fr;
    analysis_store.gamma_rr 			= gamma_rr;
    analysis_store.gamma_rl 			= gamma_rl;
    analysis_store.gamma_fl 			= gamma_fl;
    analysis_store.gamma_fr 			= gamma_fr;
    analysis_store.DFx_f 				= DFx_f;
	analysis_store.DFx_r 				= DFx_r;
    analysis_store.DFy_f 				= DFy_f;
    analysis_store.DFy_r 				= DFy_r;
	analysis_store.DFz_f 				= DFz_f;
	analysis_store.DFz_r 				= DFz_r;
    analysis_store.mu_r 				= mu_r;
    analysis_store.mu_f 				= mu_f;
    analysis_store.Ay_norm 				= Ay_norm;
    analysis_store.handling 			= handling;
    analysis_store.ay_fit_lin 			= ay_fit_lin;
    analysis_store.handling_fit_lin 	= handling_fit_lin;
    analysis_store.ay_fit_nonlin 		= ay_fit_nonlin;
    analysis_store.handling_fit_nonlin 	= handling_fit_nonlin;
    analysis_store.K_US_theo 			= K_US_theo;
    analysis_store.K_US_theo2 			= K_US_theo2;
    analysis_store.tau_D 				= tau_D;
    analysis_store.Fz_r0 				= Fz_r0;
    analysis_store.Fz_f0 				= Fz_f0;
    
end
    