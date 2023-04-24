%% FITTING TYRE DATA-LATERAL BEHAVIOUR WITH MAGIC FORMULA 96
%% Syntax functions
% plot_fitted_data(x_raw, y_raw, x_fit, y_fit, label_x, label_y, name, plot_title, line_width, font_size_title)
% plot_selected_data(tab, font_size_title)
% magic_formula_stiffness(x, B, C, D, E, SV)
% initialise_ty_data(R0, Fz0)


%% Initialisation
clc
clearvars 
close all 

addpath('..\utilities\')
addpath('..\dataset\')
addpath('\MAINS\..')
addpath('tyre_lib\FY')

% Suppress warning export_fig for docked figures
warning('off', 'MATLAB:Figure:SetPosition')

% Choice of the dataset
data_set_path = 'dataset/';
data_set = 'Goodyear_B1464run13.mat';
struct_name = 'Goodyear_B1464';  
load_type = 'lateral'; 

initialization


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

%% -------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
[TData0, ~] = intersect_table_data(GAMMA_0, FZ_220 );

% figure 3
plot_selected_data(TData0, font_size_title);
% if exist(['tyre_', struct_name,'.mat'], 'file')
%   load(['tyre_', struct_name,'.mat']);
% else
  tyre_coeffs = initialise_tyre_data(R0, Fz0);
% end

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
P0 = [  1,   2,   1,     0,   1,  0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [pCy1 pDy1 pEy1  pHy1  pKy1  pVy1]
lb = [1,  0.1,   -1,   -10,  -4, -10]; % lower bound
ub = [2,    4,   1,    10,   100, 10]; % upper bound

ALPHA_vec = TData0.SA;  % side slip angle
FY_vec    = TData0.FY;  % lateral force

% check guess
SA_vec = -0.3:0.001:0.3; % lateral slip to be used in the plot for check the interpolation outside the input data range

% resid_pure_Fy returns the residual, so minimize the residual varying Y. 
% It is an unconstrained minimization problem 
[P_fz_nom, res_Fy0, exitflag] = fmincon(@(P)resid_pure_Fy(P, FY_vec, ALPHA_vec, 0, FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);
  
%%
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

res_Fy0 = resid_pure_Fy(P_fz_nom, FY_vec, ALPHA_vec, 0, FZ0, tyre_coeffs);


%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
% extract data with variable load
[TDataDFz, ~] = GAMMA_0; %intersect_table_data(GAMMA_0);
TDataDFz([9808:9930], :) = [];


zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));

% Guess values for parameters to be optimised
%    [ pHy2, pDy2, pEy2, pVy2, pKy2, pEy3]
%% WORK HERE
P0 = 0*ones(6, 1);
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [pHy2, pDy2, pEy2, pVy2, pKy2, pEy3]
lb = -inf*ones(6, 1);
ub = inf*ones(6, 1);

ALPHA_vec = TDataDFz.SA;
FY_vec    = TDataDFz.FY;
FZ_vec    = TDataDFz.FZ;

% check guess
SA_vec = -0.3:0.001:0.3;

% LSM_pure_Fy returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_dfz, res_Fy0_varFz, exitflag] = fmincon(@(P)resid_pure_Fy_varFz(P, FY_vec, ALPHA_vec, 0, FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);
P_dfz
res_Fy0_varFz

%%
% Change tyre data with new optimal values     

tyre_coeffs.pHy2 = P_dfz(1);
tyre_coeffs.pDy2 = P_dfz(2);
tyre_coeffs.pEy2 = P_dfz(3);
tyre_coeffs.pVy2 = P_dfz(4);
tyre_coeffs.pKy2 = P_dfz(5);
tyre_coeffs.pEy3 = P_dfz(6);

res_Fy0_varFz = resid_pure_Fy_varFz(P_dfz, FY_vec, SA_vec, 0, FZ_vec, tyre_coeffs);

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));

FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1550.FZ)*tmp_ones,tyre_coeffs);

x_fit = repmat(SA_vec', 1, 5);
y_fit = [FY0_fz_var_vec1', FY0_fz_var_vec2', FY0_fz_var_vec3', FY0_fz_var_vec4', FY0_fz_var_vec5'];
data_label = string(5);
data_label = ["220 [N]", "440 [N]", "700 [N]", "1120 [N]", "1550 [N]"];
%figure 6
plot_fitted_data(TDataDFz.SA, TDataDFz.FY, x_fit, y_fit, '$\alpha$ [-]', '$F_{Y}$ [N]', data_label,'fig_fit_variable_load_FY', 'Fitting with variable load', line_width, font_size_title);


[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_440.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_1550.FZ), tyre_coeffs);
Calfa_vec5_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);

Calfa_vec1 = MF96_CorneringStiffness_FY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness_FY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness_FY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness_FY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec5 = MF96_CorneringStiffness_FY(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1550.FZ)*tmp_ones,tyre_coeffs);

% figure 7
plot_stiffness_FY; 

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER
%--------------------------------------------------------------------------
% extract data with variable load
[TDataGamma, ~] = FZ_220; % intersect_table_data(SA_0, FZ_220);

% Fit the coeffs { pDy3 pEy4 pHy3 pKy3 pVy3 pVy4}

% Guess values for parameters to be optimised
%    [pHy3, pDy3, pKy3, pEy4, pVy3, pVy4]
P0 = 0*ones(6, 1); %[ -10 10 0 0 0   ]; 

% NOTE: many local minima => limits on parameters are fundamentals
lb = -inf*ones(6, 1);
ub = inf*ones(6, 1);

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
export_fig(fig_camber_FY, 'images\fig_camber_FY.png')

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma, res_Fy0_varGamma, exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub)

%%
% Change tyre data with new optimal values
tyre_coeffs.pHy3 = P_varGamma(1);
tyre_coeffs.pDy3 = P_varGamma(2);
tyre_coeffs.pKy3 = P_varGamma(3);
tyre_coeffs.pEy4 = P_varGamma(4);
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

plot_y = cell(5,1);
plot_x = cell(5,1);
plot_fit = cell(5,1);
plot_x = {TDataGamma0.SA, TDataGamma1.SA, TDataGamma2.SA, TDataGamma3.SA, TDataGamma4.SA};
plot_y = {TDataGamma0.FY, TDataGamma1.FY, TDataGamma2.FY, TDataGamma3.FY, TDataGamma4.FY};
plot_fit = {FY0_Gamma0, FY0_Gamma1, FY0_Gamma2, FY0_Gamma3, FY0_Gamma4};
plot_label = ["$\gamma = 0 [deg]$", "$\gamma = 1 [deg]$", "$\gamma = 2 [deg]$", "$\gamma = 3 [deg]$", "$\gamma = 4 [deg]$"];

plot_fitted_data_struct(plot_x, plot_y, plot_x, plot_fit, ...
  '$\alpha$ [-]', '$F_{Y}$ [N]', plot_label, 'fig_fit_variable_camber_FY', ...
  'Fitting with variable camber', line_width, font_size_title, colors_vect);


%% -------------------
% PLOT THE RESIDUALS
fprintf("Pure confitions: %6.3f\n", res_Fy0);
fprintf("Variable load: %6.3f\n", res_Fy0_varFz);
fprintf("Variable camber: %6.3f\n", res_Fy0_varGamma);


%%
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




