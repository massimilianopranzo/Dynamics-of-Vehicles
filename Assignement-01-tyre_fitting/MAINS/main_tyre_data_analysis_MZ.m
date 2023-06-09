%% FITTING TYRE DATA-SELF ALIGNING MOMENT WITH MAGIC FORMULA 96

%% Initialisation
clc
clearvars 
close all 

addpath('utilities\')

% Choice of the dataset
data_set_path = 'dataset/';
data_set = 'Goodyear_B1464run13.mat';
struct_name = 'Goodyear_B1464';  
load_type = 'self_aligning'; 

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
cut_end = length(MZ);
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

%%-------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
if exist(['tyre_' struct_name,'.mat'], 'file')
  load(['tyre_' struct_name,'.mat']);
else
  tyre_coeffs = initialise_tyre_data(R0, Fz0);
end
tyre_coeffs.LKY = 1;
tyre_coeffs.LT = 1;
tyre_coeffs.LMR = 1;

%%
FY_vec = TData0.FY;
FZ0 = mean(TData0.FZ);
zeros_vec = zeros(size(TData0.SA));
ones_vec  = ones(size(TData0.SA));
FY0_guess = MF96_FY0_vec(zeros_vec, TData0.SA, zeros_vec, tyre_coeffs.FZ0 * ones_vec, tyre_coeffs);
MZ0_guess = MF96_MZ0_vec(TData0.SA, zeros_vec, tyre_coeffs.FZ0 * ones_vec, FY0_guess, tyre_coeffs); 
label = "Fitted data";
% figure 4
plot_fitted_data(TData0.SA, TData0.MZ, TData0.SA, MZ0_guess, '$\alpha [deg]$', '$M_{Z0}$ [N]', label, 'fig_fit_guess_MZ', 'Fitting with initial guess', line_width, font_size_title)


%% ------------------------------------------------------------------------
% FITTING WITH NOMINAL LOAD Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10
%--------------------------------------------------------------------------
% Guess values for parameters to be optimised
%    [qHz1, qBz1, qCz1, qDz1, qEz1, qEz4, qBz9, qDz6, qBz10]
P0 = [   0,0,0,0,0,0,0,0,0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [qHz1, qBz1, qCz1, qDz1, qEz1, qEz4, qBz9, qDz6, qBz10]
lb = [0, -5, -5, 0, -5, 0, 0, 0, 0, 0, -5]; % lower bound
ub = [15, 10, 10 , 10, 10, 10, 10, 50, 10 ]; % upper bound
% lb = [];
% ub = [];

ALPHA_vec = TData0.SA;  % lateral slip angle
MZ_vec    = TData0.MZ;  % aligning moment
FY_vec    = TData0.FY;  % laterla force

% check guess
SA_vec = -0.3:0.001:0.3; % lateral slip to be used in the plot for check the interpolation outside the input data range

% resid_pure_Mz returns the residual, so minimize the residual varying Y. 
% It is an unconstrained minimization problem 
[P_fz_nom, fval, exitflag] = fmincon(@(P)resid_pure_Mz(P, MZ_vec, ALPHA_vec, 0, FZ0, FY_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

P_fz_nom
fval

% Update tyre data with new optimal values                   
tyre_coeffs.qHz1  = P_fz_nom(1); 
tyre_coeffs.qBz1  = P_fz_nom(2);
tyre_coeffs.qCz1  = P_fz_nom(3);
tyre_coeffs.qDz1  = P_fz_nom(4);
tyre_coeffs.qEz1  = P_fz_nom(5);
tyre_coeffs.qEz4  = P_fz_nom(6);
tyre_coeffs.qBz9  = P_fz_nom(7);
tyre_coeffs.qDz6 = P_fz_nom(8);
tyre_coeffs.qBz10  = P_fz_nom(9);

FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), tyre_coeffs.FZ0.*ones(size(SA_vec)), tyre_coeffs);
MZ0_fz_nom_vec = MF96_MZ0_vec(SA_vec, zeros(size(SA_vec)), FZ0.*ones(size(SA_vec)), FY0_fz_nom_vec, tyre_coeffs);
        
data_label = "Fitted data"; 
% figure 5          
plot_fitted_data(TData0.SA, TData0.MZ, SA_vec', MZ0_fz_nom_vec', '$\alpha [-]$', '$M_{z0}$ [N]', data_label, 'fig_fit_pure_conditions_MZ', 'Fitting in pure conditions', line_width, font_size_title)


%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
% extract data with variable load
[TDataDFz, ~] = GAMMA_0; % intersect_table_data(SA_0, GAMMA_0);

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
%FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));

% Guess values for parameters to be optimised
%    [  qHz2 qHz3 qHz4 qBz2 qBz3 qEz2 qEz3 qDz7 qDz8 qDz9]
%% WORK HERE
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [ qHz2, qBz2, qBz3, qDz2, qEz2, qEz3, qDz7]
P0 = [0, 0, 0, 0, 0, 0, 0];
%% WORK HERE
lb = [0, 0, 0, 0, 0, 0, 0];
ub = [10, 10, 10, 10, 10, 10, 10];

ALPHA_vec = TDataDFz.SA;
MZ_vec    = TDataDFz.MZ;
FZ_vec    = TDataDFz.FZ;
FY_vec    = TDataDFz.FY;

% check guess
SA_vec = -0.3:0.001:0.3;
% FX0_dfz_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)), zeros(size(SL_vec)), ...
%                            TDataDFz.FZ,tyre_coeffs);

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Mz_varFz(P, MZ_vec, ALPHA_vec, 0, FZ_vec, FY_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values          
tyre_coeffs.qHz2 = P(1); 
tyre_coeffs.qBz2 = P(2); 
tyre_coeffs.qBz3 = P(3); 
tyre_coeffs.qDz2 = P(4); 
tyre_coeffs.qEz2 = P(5); 
tyre_coeffs.qEz3 = P(6); 
tyre_coeffs.qDz7 = P(7); 

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));
%%
FY0_fz_var_vec1 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec1 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones, FY0_fz_var_vec1, tyre_coeffs);
FY0_fz_var_vec2 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec2 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones, FY0_fz_var_vec2, tyre_coeffs);
FY0_fz_var_vec3 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec3 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones, FY0_fz_var_vec3, tyre_coeffs);
FY0_fz_var_vec4 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec4 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones, FY0_fz_var_vec4, tyre_coeffs);

x_fit = repmat(SA_vec', 1, 4);
y_fit = [MZ0_fz_var_vec1', MZ0_fz_var_vec2', MZ0_fz_var_vec3', MZ0_fz_var_vec4'];
data_label = string(4);
data_label = ["220 [N]", "700 [N]", "900 [N]", "1120 [N]"];
%figure 6
plot_fitted_data(TDataDFz.SA, TDataDFz.MZ, x_fit, y_fit, '$\alpha$ [-]', '$M_{Z}$ [N]', data_label,'fig_fit_variable_load_MZ', 'Fitting with variable load', line_width, font_size_title)

%% Stiffness
[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_220.FZ), tyre_coeffs);
% Calfa_vec1_0 = magic_formula_stiffness_MZ(alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);
[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_700.FZ), tyre_coeffs);
% Calfa_vec2_0 = magic_formula_stiffness_MZ(alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);
[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_900.FZ), tyre_coeffs);
% Calfa_vec3_0 = magic_formula_stiffness_MZ(alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);
[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_1120.FZ), tyre_coeffs);
% Calfa_vec4_0 = magic_formula_stiffness_MZ(alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);

Calfa_vec1 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_220.FZ)  * tmp_ones, FY0_fz_var_vec1, tyre_coeffs);
Calfa_vec1 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_440.FZ)  * tmp_ones, FY0_fz_var_vec2, tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_700.FZ)  * tmp_ones, FY0_fz_var_vec3, tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_1120.FZ) * tmp_ones, FY0_fz_var_vec4, tyre_coeffs);

% figure 7
% plot_stiffness_MZ; 

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER
%--------------------------------------------------------------------------
% extract data with variable load
[TDataGamma, ~] = FZ_220; % intersect_table_data(SA_0, FZ_220);

% Guess values for parameters to be optimised
%    [qHz3, qBz4, qBz5, qDz3, qDz4, qEz5, qDz8]
P0 = [1, 1, 1, 1, 1, 1, 1]; 

% NOTE: many local minima => limits on parameters are fundamentals
lb = [0, 0, 0, 0, 0, 0, 0];
ub = [10, 10, 10, 10, 10, 10, 10];

zeros_vec = zeros(size(TDataGamma.SA));
ones_vec  = ones(size(TDataGamma.SA));

ALPHA_vec = TDataGamma.SA;
GAMMA_vec = TDataGamma.IA; 
MZ_vec    = TDataGamma.MZ;
FZ_vec    = TDataGamma.FZ;
FY_vec    = TDataGamma.FY;

TDataGamma0 = intersect_table_data(GAMMA_0, FZ_220);
TDataGamma1 = intersect_table_data(GAMMA_1, FZ_220);
TDataGamma2 = intersect_table_data(GAMMA_2, FZ_220);
TDataGamma3 = intersect_table_data(GAMMA_3, FZ_220);
TDataGamma4 = intersect_table_data(GAMMA_4, FZ_220);

FY0_Gamma0 = MF96_FY0_vec(zeros_vec, TDataGamma0.SA, TDataGamma0.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma0 = MF96_MZ0_vec(TDataGamma0.SA, TDataGamma0.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma0, tyre_coeffs);

FY0_Gamma1 = MF96_FY0_vec(zeros_vec, TDataGamma1.SA, TDataGamma1.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma1 = MF96_MZ0_vec(TDataGamma1.SA, TDataGamma1.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma1, tyre_coeffs);

FY0_Gamma2 = MF96_FY0_vec(zeros_vec, TDataGamma2.SA, TDataGamma2.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma2 = MF96_MZ0_vec(TDataGamma2.SA, TDataGamma2.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma2, tyre_coeffs);

FY0_Gamma3 = MF96_FY0_vec(zeros_vec, TDataGamma3.SA, TDataGamma3.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma3 = MF96_MZ0_vec(TDataGamma3.SA, TDataGamma3.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma3, tyre_coeffs);

FY0_Gamma4 = MF96_FY0_vec(zeros_vec, TDataGamma4.SA, TDataGamma4.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma4 = MF96_MZ0_vec(TDataGamma4.SA, TDataGamma4.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma4, tyre_coeffs);

fig_fit_variable_camber_MZ = figure('Color', 'w');
hold on
plot(TDataGamma0.SA, TDataGamma0.MZ, '.', 'Color', colors_vect(1, :), 'HandleVisibility', 'off');
plot(TDataGamma0.SA, MZ0_Gamma0, '-', 'Color', colors_vect(1, :) , 'LineWidth', line_width, 'DisplayName', '$\gamma = 0 [deg]$' );

plot(TDataGamma1.SA, TDataGamma1.MZ, '.', 'Color', colors_vect(2, :), 'HandleVisibility', 'off');
plot(TDataGamma1.SA, MZ0_Gamma1, '-', 'Color', colors_vect(2, :), 'LineWidth', line_width , 'DisplayName', '$\gamma = 1 [deg]$');

plot(TDataGamma2.SA, TDataGamma2.MZ, '.', 'Color', colors_vect(3, :), 'HandleVisibility', 'off');
plot(TDataGamma2.SA, MZ0_Gamma2, '-', 'Color', colors_vect(3, :) , 'LineWidth', line_width,  'DisplayName', '$\gamma = 2 [deg]$' );

plot(TDataGamma3.SA, TDataGamma3.MZ, '.', 'Color', colors_vect(4, :), 'HandleVisibility', 'off');
plot(TDataGamma3.SA, MZ0_Gamma3, '-', 'Color', colors_vect(4, :), 'LineWidth', line_width, 'DisplayName', '$\gamma = 3 [deg]$' );

plot(TDataGamma4.SA, TDataGamma4.MZ, '.', 'Color', colors_vect(5, :), 'HandleVisibility', 'off');
plot(TDataGamma4.SA, MZ0_Gamma4, '-', 'Color', colors_vect(5, :)  , 'LineWidth', line_width, 'DisplayName', '$\gamma = 4 [deg]$' );
hold off
legend('location', 'northeast')
xlabel('$\alpha$')
ylabel('$M_{Z}$ [N]')
title('Fitting with variable camber $F_{Z}$ = 220[N]', 'FontSize',font_size_title)
grid on
export_fig(fig_fit_variable_camber_MZ, 'images\fig_fit_variable_camber_MZ.png');

% fig_camber_MZ = figure('Color', 'w');
% plot(ALPHA_vec, MZ_vec, '.', 'DisplayName', '$F_{z}=220 [N]$');
% grid on
% xlabel('$\gamma$ [rad]')
% ylabel('$M_Z$ [N]')
% title('Self aligning moment as function of camber angle')
% legend('Location', 'best')
% export_fig(fig_camber_MZ, 'images\fig_camber_MZ.png')

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma, fval, exitflag] = fmincon(@(P)resid_pure_Mz_varGamma(P, MZ_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0, FY_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values
tyre_coeffs.qBz4 = P_varGamma(1); 
tyre_coeffs.qBz5 = P_varGamma(2); 
tyre_coeffs.qEz5 = P_varGamma(3); 
%%
FY0_varGamma_vec = MF96_FY0_vec(zeros_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_varGamma_vec = MF96_MZ0_vec(ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec, FY_varGamma, tyre_coeffs);

data_label = ["Fitted data FZ=220 [N]"];
x_fit = repmat(ALPHA_vec, 1, 2);
y_fit = [MZ0_varGamma_vec];
plot_fitted_data(TDataGamma.SA, TDataGamma.MZ, x_fit, y_fit, '$\alpha$ [-]', '$M_{Z}$ [N]', data_label, 'fig_fit_variable_camber_MZ', 'Fitting with variable camber', line_width, font_size_title)

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n', 1 - res_Mz0_varGamma);

% [~, ~, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(alpha, phi, Fz, tyre_data)
% fprintf('By      = %6.3f\n', By);
% fprintf('Cy      = %6.3f\n', Cy);
% fprintf('muy     = %6.3f\n', Dy/tyre_coeffs.FZ0);
% fprintf('Ey      = %6.3f\n', Ey);
% fprintf('SVy     = %6.3f\n', SVy);
% fprintf('alpha_y = %6.3f\n', alpha__y);
% fprintf('Ky      = %6.3f\n', By*Cy*Dy/tyre_coeffs.FZ0);


%% Save tyre data structure to mat file
%
save(['tyre_' struct_name,'.mat'],'tyre_coeffs');



