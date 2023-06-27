 droda = -diff(rho_ss)./diff(Ay)*g*L; % diffentiation of the turning radious
 fig_KUS = figure('Name','Prova','NumberTitle','off', 'Color','w'); clf
  set(gca, 'FontSize',45);
  hold on
  grid on; box on;
  if sim_options.test_type == 1
    plot(Ay_hand(1:end-1), K_US_theo, 'LineWidth',2, 'DisplayName','$-\frac{mg}{L}(\frac{L_f}{K_{yr}} - \frac{Lr}{K_{yf}})$')
    plot(Ay_hand(1:end-1), K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$- \frac{d\Delta \alpha}{d a_y / g}$')
    plot(ay_fit_lin, K_US*ones(length(ay_fit_lin), 1), '--r', 'LineWidth',2, 'DisplayName', 'Linear $K_{US}$');
  elseif sim_options.test_type == 2
    plot(Ay_hand(1:end-1), K_US_theo, 'LineWidth',2, 'DisplayName','$-\frac{mg}{L^2}(\frac{L_f}{K_{yr}} - \frac{Lr}{K_{yf}})$')
    plot(Ay_hand(1:end-1), K_US_theo2, '--', 'LineWidth',2, 'DisplayName','$- \frac{1}{L}\frac{d\Delta \alpha}{d a_y / g}$')
    plot(ay_fit_lin, K_US*ones(length(ay_fit_lin), 1)/L, '--r', 'LineWidth',2, 'DisplayName', 'Linear $K_{US}$');
  end
  %plot(Ay_norm(1:end-1), K_US_theo3, 'LineWidth',2, 'DisplayName','Formula 2')
  plot(Ay(20000:end-1)/g, droda(20000:end), 'g')
  title('Normalized understeering gradient')
%   legend('location', 'southwest')
  xlim("padded")
  ylabel('g$K_{US}$')
  xlabel('ay/g [-]')
  set(gca, 'FontSize',45); 
  if enable_export == 1;
    export_figure(fig_KUS, strcat('\fig_KUS', suffix, '.eps'), 'images\');
  end