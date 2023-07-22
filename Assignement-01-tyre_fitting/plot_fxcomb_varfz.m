fig = figure('Color', 'w'); clf;
markers = ['.', 'o', '+'];
loads = ["FZ = 220 [N]","FZ = 700 [N]", "FZ = 900 [N]", "FZ = 1120 [N]"];
sgtitle('Fitting with variable load $F_{Z}$ = 220[N]', 'FontSize',font_size_title)
for i = 1:4 % {i}{j} i = load, j = alpha
    subplot(2,2,i); hold on; grid on; box on;
    if i == 4
        hvis = 'on';
    else
        hvis = 'off';
    end
    for j = 1:3
        plot(x_raw{i}{j}(:), FX_raw{i}{j}(:), markers(j),'Color', colors_vect(j, :), 'HandleVisibility', 'off', 'MarkerSize',5);
        plot(kappa_var, FX_fit{i}{j}(:), '-', 'Color', colors_vect(j, :),'HandleVisibility',hvis ,'LineWidth', 2,  'DisplayName', ['$\alpha =', num2str(2*j-2), ' [deg]$'] );
    end
    if i == 4
        legend('location', 'southeast')
    end
    title(loads(i))
    hold off
end

xlabel('$\kappa$')
ylabel('$F_{x0}$ [N]')


% 
% plot(x_raw, TDataGamma0.FX, '.', 'Color', colors_vect(1, :), 'HandleVisibility', 'off', 'MarkerSize',10);
% plot(TDataGamma0.SL, FX0_Gamma0, '-', 'Color', colors_vect(1, :) , 'LineWidth', line_width, 'DisplayName', '$\gamma = 0 [deg]$' );
%     
% plot(TDataGamma2.SL, TDataGamma2.FX, '.', 'Color', colors_vect(3, :), 'HandleVisibility', 'off', 'MarkerSize',10);
% plot(TDataGamma2.SL, FX0_Gamma2, '-', 'Color', colors_vect(3, :) , 'LineWidth', line_width,  'DisplayName', '$\gamma = 2 [deg]$' );
%    
% plot(TDataGamma4.SL, TDataGamma4.FX, '.', 'Color', colors_vect(5, :), 'HandleVisibility', 'off', 'MarkerSize',10);
% plot(TDataGamma4.SL, FX0_Gamma4, '-', 'Color', colors_vect(5, :)  , 'LineWidth', line_width, 'DisplayName', '$\gamma = 4 [deg]$' );
% 
% hold off
% legend('location', 'southeast')
% xlabel('$\kappa$')
% ylabel('$F_{x0}$ [N]')
% title('Fitting with variable camber $F_{Z}$ = 220[N]', 'FontSize',font_size_title)
% export_fig(fig, 'images\fig_fit_variable_camber_FX_1plot')