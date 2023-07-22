if strcmp(act, "longitudinal")
    fig = figure('Color', 'w');
    hold on; grid on; box on;
    plot(TDataGamma0.SL, TDataGamma0.FX, '.', 'Color', colors_vect(1, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(TDataGamma0.SL, FX0_Gamma0, '-', 'Color', colors_vect(1, :) , 'LineWidth', line_width, 'DisplayName', '$\gamma = 0 [deg]$' );
        
    plot(TDataGamma2.SL, TDataGamma2.FX, '.', 'Color', colors_vect(3, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(TDataGamma2.SL, FX0_Gamma2, '-', 'Color', colors_vect(3, :) , 'LineWidth', line_width,  'DisplayName', '$\gamma = 2 [deg]$' );
       
    plot(TDataGamma4.SL, TDataGamma4.FX, '.', 'Color', colors_vect(5, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(TDataGamma4.SL, FX0_Gamma4, '-', 'Color', colors_vect(5, :)  , 'LineWidth', line_width, 'DisplayName', '$\gamma = 4 [deg]$' );
   
    hold off
    legend('location', 'southeast')
    xlabel('$\kappa$')
    ylabel('$F_{x0}$ [N]')
    title('Fitting with variable camber $F_{Z}$ = 220[N]', 'FontSize',font_size_title)
    export_fig(fig, 'images\fig_fit_variable_camber_FX_1plot')
elseif strcmp(act, "lat")
    fig = figure('Color', 'w');
    hold on; grid on; box on;
    plot(TDataGamma0.SA, TDataGamma0.FY, '.', 'Color', colors_vect(1, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, FY0_Gamma0, '-', 'Color', colors_vect(1, :) , 'LineWidth', line_width, 'DisplayName', '$\gamma = 0 [deg]$' );
    
    plot(TDataGamma1.SA, TDataGamma1.FY, '.', 'Color', colors_vect(2, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, FY0_Gamma1, '-', 'Color', colors_vect(2, :), 'LineWidth', line_width , 'DisplayName', '$\gamma = 1 [deg]$');
    
    plot(TDataGamma2.SA, TDataGamma2.FY, '.', 'Color', colors_vect(3, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, FY0_Gamma2, '-', 'Color', colors_vect(3, :) , 'LineWidth', line_width,  'DisplayName', '$\gamma = 2 [deg]$' );
    
    plot(TDataGamma3.SA, TDataGamma3.FY, '.', 'Color', colors_vect(4, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, FY0_Gamma3, '-', 'Color', colors_vect(4, :), 'LineWidth', line_width, 'DisplayName', '$\gamma = 3 [deg]$' );
    
    plot(TDataGamma4.SA, TDataGamma4.FY, '.', 'Color', colors_vect(5, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, FY0_Gamma4, '-', 'Color', colors_vect(5, :)  , 'LineWidth', line_width, 'DisplayName', '$\gamma = 4 [deg]$' );
    hold off
    legend('location', 'northeast')
    xlabel('$\alpha$')
    ylabel('$F_{y0}$ [N]')
    title('Fitting with variable camber $F_{Z}$ = 220[N]', 'FontSize',font_size_title)
    export_fig(fig, 'images\fig_fit_variable_camber_FY_1plot')
elseif strcmp(act, "torque")
    fig = figure('Color', 'w');
    hold on; grid on; box on;
    plot(TDataGamma0.SA, TDataGamma0.MZ, '.', 'Color', colors_vect(1, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, MZ0_Gamma0, '-', 'Color', colors_vect(1, :) , 'LineWidth', line_width, 'DisplayName', '$\gamma = 0 [deg]$' );
    
    plot(TDataGamma1.SA, TDataGamma1.MZ, '.', 'Color', colors_vect(2, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, MZ0_Gamma1, '-', 'Color', colors_vect(2, :), 'LineWidth', line_width , 'DisplayName', '$\gamma = 1 [deg]$');
    
    plot(TDataGamma2.SA, TDataGamma2.MZ, '.', 'Color', colors_vect(3, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, MZ0_Gamma2, '-', 'Color', colors_vect(3, :) , 'LineWidth', line_width,  'DisplayName', '$\gamma = 2 [deg]$' );
    
    plot(TDataGamma3.SA, TDataGamma3.MZ, '.', 'Color', colors_vect(4, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, MZ0_Gamma3, '-', 'Color', colors_vect(4, :), 'LineWidth', line_width, 'DisplayName', '$\gamma = 3 [deg]$' );
    
    plot(TDataGamma4.SA, TDataGamma4.MZ, '.', 'Color', colors_vect(5, :), 'HandleVisibility', 'off', 'MarkerSize',10);
    plot(alpha_vec, MZ0_Gamma4, '-', 'Color', colors_vect(5, :)  , 'LineWidth', line_width, 'DisplayName', '$\gamma = 4 [deg]$' );
    hold off
    legend('location', 'northeast')
    xlabel('$\alpha$')
    ylabel('$M_{z0}$ [N]')
    title('Fitting with variable camber $F_{Z}$ = 220[N]', 'FontSize',font_size_title)
    export_fig(fig, 'images\fig_fit_variable_camber_MZ_1plot')
end

