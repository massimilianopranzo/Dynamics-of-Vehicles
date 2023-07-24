%% PLOT CORNERING STIFFNESS

fig_stiffness_MZ = figure('Name','C_alpha', 'Color', 'w');
subplot(2,1,1)
hold on
plot(mean(FZ_220.FZ),Calfa_vec1_0,'+','LineWidth',2, 'Color', colors_vect(2,:), 'MarkerSize',10)
plot(mean(FZ_440.FZ),Calfa_vec2_0,'+','LineWidth',2, 'Color', colors_vect(3,:), 'MarkerSize',10)
plot(mean(FZ_700.FZ),Calfa_vec3_0,'+','LineWidth',2, 'Color', colors_vect(4,:), 'MarkerSize',10)
plot(mean(FZ_900.FZ),Calfa_vec4_0,'+','LineWidth',2, 'Color', colors_vect(5,:), 'MarkerSize',10)
plot(mean(FZ_1120.FZ),Calfa_vec5_0,'+','LineWidth',2, 'Color', colors_vect(6,:), 'MarkerSize',10)
legend({'$Fz_{220}$', '$Fz_{440}$', '$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'}, 'Location', 'best')
xlabel('Vertical load [N]')
ylabel('$C_{\alpha, MAX}$ [Nm]')
grid on
title('Self aligning moment stiffness for $\alpha = 0$')

subplot(2,1,2)
hold on
plot(SA_vec,Calfa_vec1,'-','LineWidth',2, 'Color', colors_vect(2,:))
plot(SA_vec,Calfa_vec2,'-','LineWidth',2, 'Color', colors_vect(3,:))
plot(SA_vec,Calfa_vec3,'-','LineWidth',2, 'Color', colors_vect(4,:))
plot(SA_vec,Calfa_vec4,'-','LineWidth',2, 'Color', colors_vect(5,:))
plot(SA_vec,Calfa_vec5,'-','LineWidth',2, 'Color', colors_vect(6,:))
legend({'$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'}, 'Location', 'best')
xlabel('$\alpha$ [-]')
ylabel('$C_{\alpha}$ [Nm]')
grid on
title('Self aligning moment stiffness')

sgtitle('Self aligning moment stiffness', 'interpreter', 'latex', ...
  'FontSize', font_size_title)
export_fig(fig_stiffness_MZ, 'images\fig_stiffness_MZ.png')