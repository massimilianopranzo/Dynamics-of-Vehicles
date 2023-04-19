%% FITTING TYRE DATA-LATERAL BEHAVIOUR WITH MAGIC FORMULA 96


%% Initialisation
clc
clearvars 
close all 

addpath('utilities\')

% Choice of the dataset
data_set_path = 'dataset/';
data_set = 'Goodyear_B1464run13.mat';
struct_name = 'Goodyear_B1464';  

initialization

%% Syntax functions
% plot_fitted_data(x_raw, y_raw, x_fit, y_fit, label_x, label_y, name, plot_title, line_width, font_size_title)
% plot_selected_data(tab, font_size_title)
% magic_formula_stiffness(x, B, C, D, E, SV)
% initialise_ty_data(R0, Fz0)


%% Load raw data  
load ([data_set_path, data_set]); % pure lateral

% select dataset portion
cut_start = 1;
cut_end = length(FY);
smpl_range = cut_start:cut_end;

% figure 1
plots_raw_data;

%% Sort data
sorting_data;

% figure 2
plot_sorted_data;

%% Intersect tables to obtain specific sub-datasets and plot them
[TData0, ~] = intersect_table_data(GAMMA_0, FZ_220 );

% figure 3
plot_selected_data(TData0, font_size_title);

%% -------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
if exist(['tyre_', struct_name,'.mat'], 'file')
  load(['tyre_', struct_name,'.mat']);
else
  tyre_coeffs = initialise_tyre_data(R0, Fz0);
end

FZ0 = mean(TData0.FZ);
zeros_vec = zeros(size(TData0.SA));
ones_vec  = ones(size(TData0.SA));
FY0_guess = MF96_FY0_vec(zeros_vec ,TData0.SA, zeros_vec, tyre_coeffs.FZ0 * ones_vec, tyre_coeffs); 
label = "Fitted data";
% figure 4
plot_fitted_data(TData0.SA, TData0.FY, TData0.SA, FY0_guess, '$\alpha [deg]$', '$F_{Y0}$ [N]', label, 'fig_fit_guess_FY', 'Fitting with initial guess', line_width, font_size_title)


%% ------------------------------------------------------------------------
% FITTING WITH NOMINAL LOAD Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10
%--------------------------------------------------------------------------
% Guess values for parameters to be optimised
%    [pCy1 pDy1 pEy1  pHy1  pKy1  pVy1]
P0 = [  1,   2,   1,     0,   1,   0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [pCy1 pDy1 pEy1  pHy1  pKy1  pVy1] 
lb = [1,  0.1,   0,   -10,  -50,   -10]; % lower bound
ub = [2,    4,   1,    10,   50,    10]; % upper bound
% lb = [];
% ub = [];

ALPHA_vec = TData0.SA;  % side slip angle
FY_vec    = TData0.FY;  % lateral force

% check guess
SA_vec = -0.3:0.001:0.3; % lateral slip to be used in the plot for check the interpolation outside the input data range

% resid_pure_Fy returns the residual, so minimize the residual varying Y. 
% It is an unconstrained minimization problem 
[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fy(P, FY_vec, ALPHA_vec, 0, FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDy1 = P_fz_nom(2) ;  
tyre_coeffs.pEy1 = P_fz_nom(3) ;
tyre_coeffs.pHy1 = P_fz_nom(4) ; 
tyre_coeffs.pKy1 = P_fz_nom(5) ;
tyre_coeffs.pVy1 = P_fz_nom(6) ;

FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs);
        
data_label = "Fitted data"; 
% figure 5          
plot_fitted_data(TData0.SA, TData0.FY, SA_vec', FY0_fz_nom_vec', '$\alpha [-]$', '$F_{y0}$ [N]', data_label, 'fig_fit_pure_conditions_FY', 'Fitting in pure conditions', line_width, font_size_title)


%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
% extract data with variable load
[TDataDFz, ~] = GAMMA_0; %intersect_table_data(GAMMA_0);

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
%FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));

% Guess values for parameters to be optimised
%    [   pDy2     pEy2   pHy2    pKy2   pVy2    pEy3]
%% WORK HERE
P0 = [  -0.25,   -0.25,     0,     -1,     0,   0.88].*0; 
P0 = [-0.0000   -0.2018    0.0004    4.3470   -0.0071    8.9982];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [   pDy2     pEy2   pHy2    pKy2   pVy2    pEy3]
%% WORK HERE
lb = [-5 -5 -5 -5 -5 5];
ub = [2 0 5 5 5 50];

ALPHA_vec = TDataDFz.SA;
FY_vec    = TDataDFz.FY;
FZ_vec    = TDataDFz.FZ;

% check guess
SA_vec = -0.3:0.001:0.3;

% LSM_pure_Fy returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_dfz, fval, exitflag] = fmincon(@(P)resid_pure_Fy_varFz(P, FY_vec, ALPHA_vec, 0, FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values                             
tyre_coeffs.pDy2 = P_dfz(1); % 1
tyre_coeffs.pEy2 = P_dfz(2);  
tyre_coeffs.pHy2 = P_dfz(3);
tyre_coeffs.pKy2 = P_dfz(4); 
tyre_coeffs.pVy2 = P_dfz(5);
tyre_coeffs.pEy3 = P_dfz(6); 

res_FY0_dfz_vec = resid_pure_Fy_varFz(P_dfz, FY_vec, SA_vec, 0, FZ_vec, tyre_coeffs);

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));

FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

x_fit = repmat(SA_vec', 1, 4);
y_fit = [FY0_fz_var_vec1', FY0_fz_var_vec2', FY0_fz_var_vec3', FY0_fz_var_vec4'];
data_label = string(4);
data_label = ["220 [N]", "700 [N]", "900 [N]", "1120 [N]"];
%figure 6
plot_fitted_data(TDataDFz.SA, TDataDFz.FY, x_fit, y_fit, '$\alpha$ [-]', '$F_{Y}$ [N]', data_label,'fig_fit_variable_load_FY', 'Fitting with variable load', line_width, font_size_title)


[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);

Calfa_vec1 = MF96_CorneringStiffness(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

% figure 7
plot_stiffness_FY; 

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER
%--------------------------------------------------------------------------
% extract data with variable load
[TDataGamma, ~] = FZ_220; % intersect_table_data(SA_0, FZ_220);

% Fit the coeffs { pDy3 pEy4 pHy3 pKy3 pVy3 pVy4}

% Guess values for parameters to be optimised
%    [pDy3  pEy4     pHy3    pKy3   pVy3     pVy4 ]
P0 = [5.16  -4.2     -0.03   1.3    -2.92    -2.81]; 

% NOTE: many local minima => limits on parameters are fundamentals
lb = [];
ub = [];

zeros_vec = zeros(size(TDataGamma.SA));
ones_vec  = ones(size(TDataGamma.SA));

ALPHA_vec = TDataGamma.SA;
GAMMA_vec = TDataGamma.IA; 
FY_vec    = TDataGamma.FY;
FZ_vec    = TDataGamma.FZ;

fig_camber_FY = figure('Color', 'w');
plot(ALPHA_vec,FY_vec, '.', 'DisplayName', '$F_{z}=220 [N]$');
grid on
xlabel('$\gamma$ [rad]')
ylabel('$F_Y$ [N]')
title('Lateral force as function of camber angle')
legend('Location', 'best')
export_fig(fig_camber_FY, 'images\fig_camber_FY.svg')

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values
tyre_coeffs.pDy3 = P_varGamma(1); 
tyre_coeffs.pEy4 = P_varGamma(2); 
tyre_coeffs.pHy3 = P_varGamma(3); 
tyre_coeffs.pKy3 = P_varGamma(4); 
tyre_coeffs.pVy3 = P_varGamma(5); 
tyre_coeffs.pVy4 = P_varGamma(6);


FY0_varGamma_vec = MF96_FY0_vec(zeros_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% data_label = ["Fitted data FZ=220 [N]"];
% x_fit = repmat(ALPHA_vec, 1, 2);
% y_fit = [FY0_varGamma_vec];
% plot_fitted_data(TDataGamma.SA, TDataGamma.FY, x_fit, y_fit, ...
%   '$\alpha$ [-]', '$F_{Y}$ [N]', data_label, 'fig_fit_variable_camber_FY', ...
%   'Fitting with variable camber', line_width, font_size_title)

%%
TDataGamma0 = intersect_table_data(GAMMA_0, FZ_220);
TDataGamma1 = intersect_table_data(GAMMA_1, FZ_220);
TDataGamma2 = intersect_table_data(GAMMA_2, FZ_220);
TDataGamma3 = intersect_table_data(GAMMA_3, FZ_220);
TDataGamma4 = intersect_table_data(GAMMA_4, FZ_220);
FY0_Gamma0 = MF96_FY0_vec(zeros_vec, TDataGamma0.SA, TDataGamma0.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
FY0_Gamma1 = MF96_FY0_vec(zeros_vec, TDataGamma1.SA, TDataGamma1.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
FY0_Gamma2 = MF96_FY0_vec(zeros_vec, TDataGamma2.SA, TDataGamma2.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
FY0_Gamma3 = MF96_FY0_vec(zeros_vec, TDataGamma3.SA, TDataGamma3.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
FY0_Gamma4 = MF96_FY0_vec(zeros_vec, TDataGamma4.SA, TDataGamma4.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

fig_fit_variable_camber_FY = figure('Color', 'w');
hold on
plot(TDataGamma0.SA, TDataGamma0.FY, '.', 'Color', colors_vect(1, :), 'HandleVisibility', 'off');
plot(TDataGamma0.SA, FY0_Gamma0, '-', 'Color', colors_vect(1, :) , 'LineWidth', line_width, 'DisplayName', '$\gamma = 0 [deg]$' );

plot(TDataGamma1.SA, TDataGamma1.FY, '.', 'Color', colors_vect(2, :), 'HandleVisibility', 'off');
plot(TDataGamma1.SA, FY0_Gamma1, '-', 'Color', colors_vect(2, :), 'LineWidth', line_width , 'DisplayName', '$\gamma = 1 [deg]$');

plot(TDataGamma2.SA, TDataGamma2.FY, '.', 'Color', colors_vect(3, :), 'HandleVisibility', 'off');
plot(TDataGamma2.SA, FY0_Gamma2, '-', 'Color', colors_vect(3, :) , 'LineWidth', line_width,  'DisplayName', '$\gamma = 2 [deg]$' );

plot(TDataGamma3.SA, TDataGamma3.FY, '.', 'Color', colors_vect(4, :), 'HandleVisibility', 'off');
plot(TDataGamma3.SA, FY0_Gamma3, '-', 'Color', colors_vect(4, :), 'LineWidth', line_width, 'DisplayName', '$\gamma = 3 [deg]$' );

plot(TDataGamma4.SA, TDataGamma4.FY, '.', 'Color', colors_vect(5, :), 'HandleVisibility', 'off');
plot(TDataGamma4.SA, FY0_Gamma4, '-', 'Color', colors_vect(5, :)  , 'LineWidth', line_width, 'DisplayName', '$\gamma = 4 [deg]$' );
hold off
legend('location', 'northeast')
xlabel('$\alpha$')
ylabel('$F_{Y}$ [N]')
title('Fitting with variable camber $F_{Z}$ = 220[N]', 'FontSize',font_size_title)
grid on
export_fig(fig_fit_variable_camber_FY, 'images\fig_fit_variable_camber_FY.svg');




%% Calculate the residuals with the optimal solution found above
res_Fy0_varGamma  = resid_pure_Fy_varGamma(P_varGamma,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fy0_varGamma);


[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
fprintf('By      = %6.3f\n', By);
fprintf('Cy      = %6.3f\n', Cy);
fprintf('muy     = %6.3f\n', Dy/tyre_coeffs.FZ0);
fprintf('Ey      = %6.3f\n', Ey);
fprintf('SVy     = %6.3f\n', SVy);
fprintf('alpha_y = %6.3f\n', alpha__y);
fprintf('Ky      = %6.3f\n', By*Cy*Dy/tyre_coeffs.FZ0);


%% Save tyre data structure to mat file
%
save(['tyre_' struct_name,'.mat'],'tyre_coeffs');




