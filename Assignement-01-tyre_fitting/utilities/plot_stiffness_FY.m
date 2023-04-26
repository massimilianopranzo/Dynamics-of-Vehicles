function plot_stiffness_FY(SA_vec,FZ_220, FZ_440, FZ_700, FZ_900, FZ_1120,  ...
  Calfa_vec1_0,Calfa_vec2_0, Calfa_vec3_0, Calfa_vec4_0, Calfa_vec5_0, ...
  Calfa_vec1, Calfa_vec2, Calfa_vec3, Calfa_vec4, Calfa_vec5, ...
  plot_title, plot_name, font_size_title )
%% PLOT CORNERING STIFFNESS

fig_stiffness_FY = figure('Name','C_alpha', 'Color', 'w');
subplot(2,1,1)
hold on
plot(mean(FZ_220.FZ),Calfa_vec1_0,'+','LineWidth',2, 'MarkerSize',10)
plot(mean(FZ_440.FZ),Calfa_vec2_0,'+','LineWidth',2, 'MarkerSize',10)
plot(mean(FZ_700.FZ),Calfa_vec3_0,'+','LineWidth',2, 'MarkerSize',10)
plot(mean(FZ_900.FZ),Calfa_vec4_0,'+','LineWidth',2, 'MarkerSize',10)
plot(mean(FZ_1120.FZ),Calfa_vec5_0,'+','LineWidth',2, 'MarkerSize',10)
legend({'$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'...
  }, 'Location', 'northwest')
xlabel('Vertical load [N]')
ylabel('$C_{\alpha}$ [N]')
grid on
title('Cornering stiffness for $\alpha = 0$')

subplot(2,1,2)
hold on
plot(SA_vec,Calfa_vec1,'-','LineWidth',2)
plot(SA_vec,Calfa_vec2,'-','LineWidth',2)
plot(SA_vec,Calfa_vec3,'-','LineWidth',2)
plot(SA_vec,Calfa_vec4,'-','LineWidth',2)
plot(SA_vec,Calfa_vec5,'-','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$', ...
  '$Fz_{1120}$'}, 'Location', 'northwest')
xlabel('$\alpha$ [-]')
ylabel('$C_{\alpha}$ [N]')
grid on
title('Cornering stiffness')

sgtitle(plot_title, 'Interpreter', 'latex', 'FontSize', font_size_title)
name = ['images\', plot_name, '.png'];
export_fig(fig_stiffness_FY, name)