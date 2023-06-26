const_speed = load("results_constspeed.mat");
const_steer = load("results_conststeer.mat");


fig_load_transf = figure('Name','Load transfer','NumberTitle','off', 'Color', 'w'); clf
ax(1) = subplot(221);
plot(Ay_norm, DFx_f(1:end-1),'LineWidth',2)
hold on
plot(Ay_norm, DFx_r(1:end-1), '--', 'LineWidth',2)
legend('$\Delta F_{xf}$','$\Delta F_{xr}$','location','northwest')
title('$\Delta F_{xf}$ and $\Delta F_{xr}$ [N]')
xlabel('$a_y/g [-]$')
grid on; box on;
set(gca, 'FontSize',30);
% --- DeltaFyf DeltaFyr -- %
ax(2) = subplot(222);
plot(Ay_norm, DFy_f(1:end-1),'LineWidth',2)
hold on
plot(Ay_norm, DFy_r(1:end-1), '--', 'LineWidth',2)
legend('$\Delta F_{yf}$','$\Delta F_{yr}$','location','northwest')
title('$\Delta F_{yf}$ and $\Delta F_{yr}$ [N]')
xlabel('$a_y/g [-]$')
grid on; box on;
set(gca, 'FontSize',30);
% --- DeltaFzf DeltaFzr -- %
ax(3) = subplot(223.5);
plot(Ay_norm, DFz_f(1:end-1),'LineWidth',2)
hold on
plot(Ay_norm, DFz_r(1:end-1), '--', 'LineWidth',2)
legend('$\Delta F_{zf}$','$\Delta F_{zr}$','location','northwest')
title('$\Delta F_{zf}$ and $\Delta F_{zr}$ [N]')
% xlabel('$a_y/g [-]$')
grid on; box on;
sgtitle('Lateral load transfer', 'FontSize', 25)
set(gca, 'FontSize',30); 
if enable_export == 1;
  export_figure(fig_load_transf, strcat('\fig_load_transf', suffix, '.eps'), 'images\');
end
clear ax
