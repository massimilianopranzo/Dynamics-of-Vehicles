function [] = plot_sorted_data(tyre_data, idx, vec_samples, GAMMA_0, GAMMA_1, ...
  GAMMA_2, GAMMA_3, GAMMA_4, GAMMA_5, FZ_220, FZ_440, FZ_700, FZ_900, FZ_1120, ...
  FZ_1550, load_type, SA_0, SA_3neg, SA_6neg, font_size_title, plot_title, ...
  plot_name)
%% PLOT SORTED DATA
to_deg = pi/180;

fig_sorted = figure('Color','w');

if strcmp(load_type, 'longitudinal')
  tiledlayout(3,1)
else
  tiledlayout(2,1)
end

sgtitle(plot_title,'fontsize', font_size_title, 'interpreter','latex') 
ax_list(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle', 'FontSize',font_size_title)
xlabel('Samples [-]')
ylabel('[deg]')
grid on

ax_list(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
plot(vec_samples(idx.FZ_1550),FZ_1550.FZ,'.');
title('Vertical force', 'FontSize',font_size_title)
xlabel('Samples [-]')
ylabel('[N]')
grid on

if strcmp(load_type, 'longitudinal')
  
  ax_list(3) = nexttile;
  plot(tyre_data.SA*to_deg)
  hold on
  plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
  plot(vec_samples(idx.SA_3neg),SA_3neg.SA*to_deg,'.');
  plot(vec_samples(idx.SA_6neg),SA_6neg.SA*to_deg,'.');
  title('Slide slip', 'FontSize',font_size_title)
  xlabel('Samples [-]')
  ylabel('[deg]')
  grid on
end

name = ['images/', plot_name, '.png'];
export_fig(fig_sorted, name);
