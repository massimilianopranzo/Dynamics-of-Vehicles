function plotMultiAnalysis(analysis_store, i)
  n_sim = length(analysis_store);

  f_Steer 			= figure('Name','Steer','NumberTitle','off');
  f_LatSlips 		= figure('Name','LatSlips','NumberTitle','off');
  f_LongSlips 		= figure('Name','LongSlips','NumberTitle','off');
  f_WheelCamber 	= figure('Name','WheelCamber','NumberTitle','off'); 
  f_LoadTransfer 	= figure('Name','LoadTransfer','NumberTitle','off');
  f_LongAxle 		= figure('Name','LongAxle','NumberTitle','off');
  f_LatAxle 		= figure('Name','LatAxle','NumberTitle','off');
  f_NormAxle 		= figure('Name','NormAxle','NumberTitle','off');
  f_Handling 		= figure('Name','Handling','NumberTitle','off');
  f_UnderGrad 		= figure('Name','UnderGrad','NumberTitle','off');

for i = 1:n_sim
    time_sim 			= analysis_store{i}.time_sim;
    delta_D 			= analysis_store{i}.delta_D;
    delta_fr 			= analysis_store{i}.delta_fr;
    delta_fl 			= analysis_store{i}.delta_fl;
    alpha_rr 			= analysis_store{i}.alpha_rr;
    alpha_rl 			= analysis_store{i}.alpha_rl;
    alpha_fl 			= analysis_store{i}.alpha_fl;
    alpha_fr 			= analysis_store{i}.alpha_fr;
    alpha_r 			= analysis_store{i}.alpha_r;
    alpha_f 			= analysis_store{i}.alpha_f;
    Fy_rr 				= analysis_store{i}.Fy_rr;
    Fy_rl 				= analysis_store{i}.Fy_rl;
    Fy_fl 				= analysis_store{i}.Fy_fl;
    Fy_fr 				= analysis_store{i}.Fy_fr;
    kappa_rr 			= analysis_store{i}.kappa_rr;
    kappa_rl 			= analysis_store{i}.kappa_rl;
    kappa_fl 			= analysis_store{i}.kappa_fl;
    kappa_fr 			= analysis_store{i}.kappa_fr;
    Fx_rr 				= analysis_store{i}.Fx_rr;
    Fx_rl 				= analysis_store{i}.Fx_rl;
    Fx_fl 				= analysis_store{i}.Fx_fl;
    Fx_fr 				= analysis_store{i}.Fx_fr;
    gamma_rr 			= analysis_store{i}.gamma_rr;
    gamma_rl 			= analysis_store{i}.gamma_rl;
    gamma_fl 			= analysis_store{i}.gamma_fl;
    gamma_fr 			= analysis_store{i}.gamma_fr;
    DFx_f 				= analysis_store{i}.DFx_f;
	DFx_r	 			= analysis_store{i}.DFx_r;
    DFy_f 				= analysis_store{i}.DFy_f;
    DFy_r 				= analysis_store{i}.DFy_r;
	DFz_r 				= analysis_store{i}.DFz_r;
	DFz_f 				= analysis_store{i}.DFz_f;
    mu_r 				= analysis_store{i}.mu_r;
    mu_f 				= analysis_store{i}.mu_f;
    Ay_norm 			= analysis_store{i}.Ay_norm;
    handling 			= analysis_store{i}.handling;
    ay_fit_lin 			= analysis_store{i}.ay_fit_lin;
    handling_fit_lin 	= analysis_store{i}.handling_fit_lin;
    ay_fit_nonlin 		= analysis_store{i}.ay_fit_nonlin;
    handling_fit_nonlin = analysis_store{i}.handling_fit_nonlin;
    K_US_theo 			= analysis_store{i}.K_US_theo;
    K_US_theo2 			= analysis_store{i}.K_US_theo2;
    tau_D 				= analysis_store{i}.tau_D;
    Fz_r0 				= analysis_store{i}.Fz_r0;
    Fz_f0 				= analysis_store{i}.Fz_f0;

  
  figure(f_Steer)
  % --- delta_0 --- %
  ax(1) = subplot(221);
  hold on
  plot(time_sim,delta_D,'LineWidth',2)
  grid on
  title('$\delta_0$ [deg]')
  xlim([0 time_sim(end)])
  % --- delta_fr --- %
  ax(2) = subplot(222);
  hold on
  plot(time_sim,delta_fr,'LineWidth',2)
  grid on
  title('$\delta_{fr}$ [deg]')
  xlim([0 time_sim(end)])
  % --- delta_fl --- %
  ax(3) = subplot(223);
  hold on
  plot(time_sim,delta_fl,'LineWidth',2)
  grid on
  title('$\delta_{fl}$ [deg]')
  xlim([0 time_sim(end)])
  % --- comparison --- %
  ax(4) = subplot(224);
  hold on
  plot(time_sim,delta_D/tau_D,'LineWidth',2)
  plot(time_sim,delta_fr,'LineWidth',2)
  plot(time_sim,delta_fl,'LineWidth',2)
  grid on
  legend('$\delta_D/\tau_D$','$\delta_{fr}$','$\delta_{fl}$','location','best')
  xlim([0 time_sim(end)])

  % -------------------------------
  %% Plot lateral tire slips and lateral forces
  % -------------------------------
  figure(f_LatSlips)
  % --- alpha_rr --- %
  ax(1) = subplot(331);
  hold on
  plot(time_sim,alpha_rr,'LineWidth',2)
  grid on
  title('$\alpha_{rr}$ [deg]')
  xlim([0 time_sim(end)])
  % --- alpha_rl --- %
  ax(2) = subplot(332);
  hold on
  plot(time_sim,alpha_rl,'LineWidth',2)
  grid on
  title('$\alpha_{rl}$ [deg]')
  xlim([0 time_sim(end)])
  % --- alpha_fr --- %
  ax(3) = subplot(333);
  hold on
  plot(time_sim,alpha_fr,'LineWidth',2)
  grid on
  title('$\alpha_{fr}$ [deg]')
  xlim([0 time_sim(end)])
  % --- alpha_fl --- %
  ax(4) = subplot(334);
  hold on
  plot(time_sim,alpha_fl,'LineWidth',2)
  grid on
  title('$\alpha_{fl}$ [deg]')
  xlim([0 time_sim(end)])
  % --- Fy_rr --- %
  ax(5) = subplot(335);
  hold on
  plot(time_sim,Fy_rr,'LineWidth',2)
  grid on
  title('$Fy_{rr}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fy_rl --- %
  ax(6) = subplot(336);
  hold on
  plot(time_sim,Fy_rl,'LineWidth',2)
  grid on
  title('$Fy_{rl}$ [Nm]')
  xlim([0 time_sim(end)])
  % --- Fy_fr --- %
  ax(7) = subplot(337);
  hold on
  plot(time_sim,Fy_fr,'LineWidth',2)
  grid on
  title('$Fy_{fr}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fy_fl --- %
  ax(8) = subplot(338);
  hold on
  plot(time_sim,Fy_fl,'LineWidth',2)
  grid on
  title('$Fy_{fl}$ [N]')
  xlim([0 time_sim(end)])

  % linkaxes(ax,'x')
  clear ax

  
  % ---------------------------------
  %% Plot longitudinal tire slips and longitudinal forces
  % ---------------------------------
  figure(f_LongSlips)
  % --- kappa_rr --- %
  ax(1) = subplot(331);
  hold on
  plot(time_sim,kappa_rr,'LineWidth',2)
  grid on
  title('$\kappa_{rr}$ [-]')
  xlim([0 time_sim(end)])
  % --- kappa_rl --- %
  ax(2) = subplot(332);
  hold on
  plot(time_sim,kappa_rl,'LineWidth',2)
  grid on
  title('$\kappa_{rl}$ [-]')
  xlim([0 time_sim(end)])
  % --- kappa_fr --- %
  ax(3) = subplot(333);
  hold on
  plot(time_sim,kappa_fr,'LineWidth',2)
  grid on
  title('$\kappa_{fr}$ [-]')
  xlim([0 time_sim(end)])
  % --- kappa_fl --- %
  ax(4) = subplot(334);
  hold on
  plot(time_sim,kappa_fl,'LineWidth',2)
  grid on
  title('$\kappa_{fl}$ [-]')
  xlim([0 time_sim(end)])
  % --- Fx_rr --- %
  ax(5) = subplot(335);
  hold on
  plot(time_sim,Fx_rr,'LineWidth',2)
  grid on
  title('$Fx_{rr}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fx_rl --- %
  ax(6) = subplot(336);
  hold on
  plot(time_sim,Fx_rl,'LineWidth',2)
  grid on
  title('$Fx_{rl}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fx_fr --- %
  ax(7) = subplot(337);
  hold on
  plot(time_sim,Fx_fr,'LineWidth',2)
  grid on
  title('$Fx_{fr}$ [N]')
  xlim([0 time_sim(end)])
  % --- Fx_fl --- %
  ax(8) = subplot(338);
  hold on
  plot(time_sim,Fx_fl,'LineWidth',2)
  grid on
  title('$Fx_{fl}$ [N]')
  xlim([0 time_sim(end)])
  
  % linkaxes(ax,'x')
  clear ax
  
  % ---------------------------------
  %% Plot wheel camber
  % ---------------------------------
  figure(f_WheelCamber)
  % --- gamma_rr --- %
  ax(1) = subplot(221);
  hold on
  plot(time_sim,gamma_rr,'LineWidth',2)
  grid on
  title('$\gamma_{rr}$ [deg]')
  xlim([0 time_sim(end)])
  % --- gamma_rl --- %
  ax(2) = subplot(222);
  hold on
  plot(time_sim,gamma_rl,'LineWidth',2)
  grid on
  title('$\gamma_{rl}$ [deg]')
  xlim([0 time_sim(end)])
  % --- gamma_fr --- %
  ax(3) = subplot(223);
  hold on
  plot(time_sim,gamma_fr,'LineWidth',2)
  grid on
  title('$\gamma_{fr}$ [deg]')
  xlim([0 time_sim(end)])
  % --- gamma_fl --- %
  ax(4) = subplot(224);
  hold on
  plot(time_sim,gamma_fl,'LineWidth',2)
  grid on
  title('$\gamma_{fl}$ [deg]')
  xlim([0 time_sim(end)])

  % linkaxes(ax,'x')
  clear ax
  
  % ---------------------------------
  %% Plot load transfer
  % ---------------------------------
  % --- DeltaFxf DeltaFxr -- %
  figure(f_LoadTransfer)
  ax(1) = subplot(221);
  hold on
  plot(time_sim, DFx_f,'LineWidth',2)
  plot(time_sim, DFx_r, '--', 'LineWidth',2)
  legend('$\Delta F_{xf}$','$\Delta F_{xr}$','location','best')
  title('$\Delta F_{xf}$ and $\Delta F_{xr}$ [N]')
  grid on
  % --- DeltaFyf DeltaFyr -- %
  ax(2) = subplot(222);
  hold on
  plot(time_sim, DFy_f,'LineWidth',2)
  plot(time_sim, DFy_r, '--', 'LineWidth',2)
  legend('$\Delta F_{yf}$','$\Delta F_{yr}$','location','best')
  title('$\Delta F_{yf}$ and $\Delta F_{yr}$ [N]')
  grid on
  % --- DeltaFzf DeltaFzr -- %
  ax(3) = subplot(223);
  hold on
  plot(time_sim, DFz_f,'LineWidth',2)
  plot(time_sim, DFz_r, '--', 'LineWidth',2)
  legend('$\Delta F_{zf}$','$\Delta F_{zr}$','location','best')
  title('$\Delta F_{zf}$ and $\Delta F_{zr}$ [N]')
  grid on
  sgtitle('Lateral load transfer', 'FontSize', 20)
  
  clear ax

  % ---------------------------------
  %% Plot axle characteristics
  % ---------------------------------
  figure(f_LongAxle)
  subplot(1,2,1)
  hold on
  grid on
  plot(alpha_r, Fx_rr, 'LineWidth',2)
  plot(alpha_r, Fx_rl, 'LineWidth',2)
  title('$F_{xr}$ [N]')
  xlabel('$\alpha_r$')
  ylabel('$F_{xrr}, F_{xrl}$')
  legend('$F_{xrr}$','$F_{xrl}$','location','southeast')
  xlim([0.001 0.06])

  subplot(1,2,2)
  hold on
  grid on
  plot(alpha_f, Fx_fr, 'LineWidth',2)
  plot(alpha_f, Fx_fl, 'LineWidth',2)
  title('$F_{xf}$ [N]')
  xlabel('$\alpha_f$')
  ylabel('$F_{xfr}, F_{xfl}$')
  legend('$F_{xfr}$','$F_{xfl}$','location','southeast')
  xlim([0.001 0.06])

  figure(f_LatAxle)
  subplot(1,2,1)
  hold on
  grid on
  plot(alpha_r, Fy_rr/Fz_r0, 'LineWidth',2)
  plot(alpha_r, Fy_rl/Fz_r0, 'LineWidth',2)
  title('$F_{yr}$ [N]')
  xlabel('$\alpha_r$')
  ylabel('$F_{yrr}, F_{yrl}$')
  legend('$F_{yrr}$','$F_{yrl}$','location','southeast')
  % xlim([0.001 0.06])
  subplot(1,2,2)
  hold on
  grid on
  plot(alpha_f, Fy_fr/Fz_f0, 'LineWidth',2)
  plot(alpha_f, Fy_fl/Fz_f0, 'LineWidth',2)
  title('$F_{yf}$ [N]')
  xlabel('$\alpha_f$')
  ylabel('$F_{yfr}, F_{yfl}$')
  legend('$F_{yfr}$','$F_{yfl}$','location','southeast')
  % xlim([0.001 0.06])

  % ---------------------------------
  %% Plot normalized axle characteristics
  % ---------------------------------
  % --- mu_r -- %
  figure(f_NormAxle)
  hold on
  grid on
  plot(alpha_r, mu_r, 'LineWidth',2)
  plot(alpha_f, mu_f, 'LineWidth',2)
  title('$\mu_r, \mu_f$')
  xlabel('$\alpha_r, \alpha_f$')
  ylabel('$\mu_r, \mu_f$')
  legend('$\mu_r$','$\mu_f$','location','best')
  

  % ---------------------------------
  %% Plot handling digram
  % ---------------------------------
  figure(f_Handling)
  hold on
  grid on
  plot(Ay_norm, handling, 'LineWidth',2)
  plot(ay_fit_lin, handling_fit_lin, '--', 'LineWidth',2)
  plot(ay_fit_nonlin, handling_fit_nonlin, '--', 'LineWidth',2)
  title('Handling diagram')
  ylabel('$\delta_{D}\tau_{H} - \rho L [rad]$')
  legend('Data', 'Fit in linear range', 'Fit in non-linear range', 'location', 'northeast')
  
  % ---------------------------------
  %% Plot understeering gradient
  % ---------------------------------
  figure(f_UnderGrad)
  hold on
  grid on
  plot(Ay_norm, K_US_theo, 'LineWidth',2)
  plot(Ay_norm(1:end-1), K_US_theo2, 'LineWidth',2)
  title('Understeering gradient')
  legend('Formula', 'Diff', 'location', 'northeast')
  xlim("padded")
  ylabel('$K_{US}$')

end

end