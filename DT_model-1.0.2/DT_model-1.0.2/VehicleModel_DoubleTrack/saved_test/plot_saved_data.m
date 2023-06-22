%% PLOT HANDLING AND KUS FOR DIFFERENT CONDITIONS
enable_export = 0; % 1 for save the plot, 0 for not
s0_t0_c0 = load("saved_test\s100_t0_c0.mat").data;
sp10_t0_c0 = load("saved_test\s110_t0_c0.mat").data;
sn10_t0_c0 = load("saved_test\s90_t0_c0.mat").data;
sp20_t0_c0 = load("saved_test\s120_t0_c0.mat").data;
sn20_t0_c0 = load("saved_test\s80_t0_c0.mat").data;
s0_t0_cp1 = load("saved_test\s100_t0_cp1.mat").data;
s0_t0_cp2 = load("saved_test\s100_t0_cp2.mat").data;
s0_t0_cn1 = load("saved_test\s100_t0_cn1.mat").data;
s0_t0_cn2 = load("saved_test\s100_t0_cn2.mat").data;
s0_tp1_c0 = load("saved_test\s100_tp1_c0.mat").data;
s0_tp2_c0 = load("saved_test\s100_tp2_c0.mat").data;
s0_tn1_c0 = load("saved_test\s100_tn1_c0.mat").data;
s0_tn2_c0 = load("saved_test\s100_tn2_c0.mat").data;

%% HANDLING
fig_handling_variable = figure('Color','w');
subplot(221)
hold on; grid on; box on;
plot(s0_t0_c0.Ay_hand, s0_t0_c0.handling, 'LineWidth',2, 'DisplayName','Ref.')
plot(sp10_t0_c0.Ay_hand, sp10_t0_c0.handling, '--', 'LineWidth',2, 'DisplayName','$+ 10\%$')
plot(sn10_t0_c0.Ay_hand, sn10_t0_c0.handling, '--', 'LineWidth',2, 'DisplayName','$- 10\%$')
plot(sp20_t0_c0.Ay_hand, sp20_t0_c0.handling, '--', 'LineWidth',2, 'DisplayName','$+ 20\%$')
plot(sn20_t0_c0.Ay_hand, sn20_t0_c0.handling, '--', 'LineWidth',2, 'DisplayName','$- 20\%$')
legend('location','southwest')
title('Variable front stiffness, $\delta=0[deg]$, $\gamma=0 [deg]$')
xlabel('ay/g [-]')
ylabel('$\delta_{D}\tau_{H} - \rho L \ [rad]$')
subplot(222)
hold on; grid on; box on;
plot(s0_t0_c0.Ay_hand, s0_t0_c0.handling, 'LineWidth',2, 'DisplayName','Ref')
plot(s0_t0_cp1.Ay_hand, s0_t0_cp1.handling, '--', 'LineWidth',2, 'DisplayName','$\gamma = +1 [deg]$')
plot(s0_t0_cn1.Ay_hand, s0_t0_cn1.handling, '--', 'LineWidth',2, 'DisplayName','$\gamma = -1 [deg]$')
plot(s0_t0_cp2.Ay_hand, s0_t0_cp2.handling, '--', 'LineWidth',2, 'DisplayName','$\gamma = +2 [deg]$')
plot(s0_t0_cn2.Ay_hand, s0_t0_cn2.handling, '--', 'LineWidth',2, 'DisplayName','$\gamma = -2 [deg]$')
legend('location','southwest')
title('Variable front camber, $\delta=0$ [deg]')
xlabel('ay/g [-]')
ylabel('$\delta_{D}\tau_{H} - \rho L \ [rad]$')
subplot(223)
hold on; grid on; box on;
plot(s0_t0_c0.Ay_hand, s0_t0_c0.handling, 'LineWidth',2, 'DisplayName','Ref.')
plot(s0_tp1_c0.Ay_hand, s0_tp1_c0.handling, '--', 'LineWidth',2, 'DisplayName','$\delta = +1 [deg]$')
plot(s0_tp2_c0.Ay_hand, s0_tp2_c0.handling, '--', 'LineWidth',2, 'DisplayName','$\delta = +2 [deg]$')
plot(s0_tn1_c0.Ay_hand, s0_tn1_c0.handling, '--', 'LineWidth',2, 'DisplayName','$\delta = -1 [deg]$')
plot(s0_tn2_c0.Ay_hand, s0_tn2_c0.handling, '--', 'LineWidth',2, 'DisplayName','$\delta = -2 [deg]$')
legend('location','northeast')
title('Variable front toe, $\gamma = 0$ [deg]')
xlabel('ay/g [-]')
ylabel('$\delta_{D}\tau_{H} - \rho L \ [rad]$')
if enable_export == 1
    export_figure(fig_handling_variable, '\fig_handling_variable.eps', 'images\');
end

%% UNDERSTEERING GRADIENT
fig_KUS_variable = figure('Color','w');
subplot(221);
hold on; grid on; box on;
plot(s0_t0_c0.Ay_hand(1:end-1), s0_t0_c0.K_US_theo2, 'LineWidth',2, 'DisplayName','Ref.')
plot(sp10_t0_c0.Ay_hand(1:end-1), sp10_t0_c0.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$+ 10\%$')
plot(sn10_t0_c0.Ay_hand(1:end-1), sn10_t0_c0.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$- 10\%$')
plot(sp20_t0_c0.Ay_hand(1:end-1), sp20_t0_c0.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$+ 20\%$')
plot(sn20_t0_c0.Ay_hand(1:end-1), sn20_t0_c0.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$- 20\%$')
legend('location','southwest')
title('Variable front stiffness, $\delta=0$ [deg], $\gamma=0$ [deg]')
xlabel('ay/g [-]')
ylabel('$K_{US}$')
subplot(222);
hold on; grid on; box on;
plot(s0_t0_c0.Ay_hand(1:end-1), s0_t0_c0.K_US_theo2, 'LineWidth',2, 'DisplayName','Ref.')
plot(s0_t0_cp1.Ay_hand(1:end-1), s0_t0_cp1.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$\delta=+1$ [deg]')
plot(s0_t0_cp2.Ay_hand(1:end-1), s0_t0_cp2.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$\delta=+2$ [deg]')
plot(s0_t0_cn1.Ay_hand(1:end-1), s0_t0_cn1.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$\delta=-1$ [deg]')
plot(s0_t0_cn2.Ay_hand(1:end-1), s0_t0_cn2.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$\delta=-2$ [deg]')
legend('location','southwest')
title('Variable front camber, $\delta=0$ [deg]')
xlabel('ay/g [-]')
ylabel('$K_{US}$')
subplot(223);
hold on; grid on; box on;
plot(s0_t0_c0.Ay_hand(1:end-1), s0_t0_c0.K_US_theo2, 'LineWidth',2, 'DisplayName','Ref.')
plot(s0_tp1_c0.Ay_hand(1:end-1), s0_tp1_c0.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$\gamma=+1$ [deg]')
plot(s0_tp2_c0.Ay_hand(1:end-1), s0_tp2_c0.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$\gamma=+2$ [deg]')
plot(s0_tn1_c0.Ay_hand(1:end-1), s0_tn1_c0.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$\gamma=-1$ [deg]')
plot(s0_tn2_c0.Ay_hand(1:end-1), s0_tn2_c0.K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$\gamma=-2$ [deg]')
legend('location','northeast')
title('Variable front toe, $\gamma=0$ [deg]')
xlabel('ay/g [-]')
ylabel('$K_{US}$')
if enable_export == 1
    export_figure(fig_KUS_variable, '\fig_KUS_variable.eps', 'images\');
end