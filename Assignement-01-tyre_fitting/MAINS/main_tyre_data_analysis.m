%% FITTING TYRE DATA WITH MAGIC FORMULA 96

%% Initialisation
clc
clearvars 
close all 

addpath('utilities\')
addpath('\MAINS\..')
addpath('tyre_lib\FX')

% Choice of the dataset
data_set_path = 'dataset/';
data_set = 'Hoosier_B1464run30.mat';
struct_name = 'Hoosier_B1464';  
load_type = 'longitudinal';

initialization

%% Syntax functions
% plot_fitted_data(x_raw, y_raw, x_fit, y_fit, label_x, label_y, name, plot_title, line_width, font_size_title)
% plot_selected_data(tab, font_size_title)
% magic_formula_stiffness(x, B, C, D, E, SV)
% initialise_ty_data(R0, Fz0)


%% Load raw data  
load ([data_set_path, data_set]); % pure lateral

% select dataset portion
cut_start = 19030;
cut_end = 38000;
smpl_range = cut_start:cut_end;

% figure 1
plots_raw_data(cut_start,cut_end,FZ,IA,SA,SL,P,TSTC,TSTI,TSTO, ...
  font_size_title,'raw_data_FX');

%% Sort data
sorting_data;
% figure 2
plot_sorted_data(tyre_data, idx, vec_samples, GAMMA_0, GAMMA_1, ...
  GAMMA_2, GAMMA_3, GAMMA_4, GAMMA_5, FZ_220, FZ_440, FZ_700, FZ_900, ...
  FZ_1120, FZ_1550, load_type, SA_0, SA_3neg, SA_6neg, font_size_title, ...
  'Sorted data longitudinal','sorted_data_FX');

%% Intersect tables to obtain specific sub-datasets and plot them
[TData0, ~] = intersect_table_data(SA_0, GAMMA_0, FZ_220 );

% figure 3
plot_selected_data(TData0, font_size_title);

%%-------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
if exist(['tyre_' struct_name,'.mat'], 'file')
  load(['tyre_' struct_name,'.mat']);
  tyre_coeffs.R0 = R0;
  tyre_coeffs.pHx1 = 0;
  tyre_coeffs.pCx1 = 0;
  tyre_coeffs.pDx1 = 0;
  tyre_coeffs.pKx1 = 0;
  tyre_coeffs.pEx1 = 0;
  tyre_coeffs.pEx4 = 0;
  tyre_coeffs.pVx1 = 0;
  tyre_coeffs.pHx2 = 0;
  tyre_coeffs.pDx2 = 0;
  tyre_coeffs.pKx2 = 0;
  tyre_coeffs.pKx3 = 0;
  tyre_coeffs.pEx2 = 0;
  tyre_coeffs.pEx3 = 0;
  tyre_coeffs.pVx2 = 0;
  tyre_coeffs.pDx3 = 0;
else
  tyre_coeffs = initialise_tyre_data(R0, Fz0);
end

FZ0 = mean(TData0.FZ);
zeros_vec = zeros(size(TData0.SL));
ones_vec  = ones(size(TData0.SL));

%% ------------------------------------------------------------------------
% FITTING WITH NOMINAL LOAD Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10
%--------------------------------------------------------------------------
% Guess values for parameters to be optimised
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1]
P0 = [  1,   2,   1,  0,   0,   1,   0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
lb = [1,  0.1,   0, 0, -10, 0, -10]; % lower bound
ub = [2,    4,   1,    1,  10,  100,  10]; % upper bound

KAPPA_vec = TData0.SL;  % slip ratio
FX_vec    = TData0.FX;  % longitudianl force

% check guess
SL_vec = -0.3:0.001:0.3; % longitudinal slip to be used in the plot for check the interpolation outside the input data range

% resid_pure_Fx returns the residual, so minimize the residual varying X. 
% It is an unconstrained minimization problem 
[P_fz_nom,fval_nom,exitflag] = fmincon(@(P)resid_pure_Fx(P, FX_vec, ...
  KAPPA_vec, 0, FZ0, tyre_coeffs),P0,[],[],[],[],lb,ub);
P_fz_nom
res_Fx0 = fval_nom
RMS_Fx0 = sqrt(fval_nom*sum(FX_vec.^2)/length(FX_vec));
R2_Fx0 = 1 - res_Fx0*sum(FX_vec.^2)/sum((FX_vec - mean(FX_vec)).^2);

% Update tyre data with new optimal values                             
tyre_coeffs.pCx1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDx1 = P_fz_nom(2) ;  
tyre_coeffs.pEx1 = P_fz_nom(3) ;
tyre_coeffs.pEx4 = P_fz_nom(4) ;
tyre_coeffs.pHx1 = P_fz_nom(5) ; 
tyre_coeffs.pKx1 = P_fz_nom(6) ;
tyre_coeffs.pVx1 = P_fz_nom(7) ;

FX0_fz_nom_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)), zeros(size(SL_vec)), ...
                              FZ0.*ones(size(SL_vec)),tyre_coeffs);
        
% figure 5   
data_label = "Fitted data"; 
plot_fitted_data(TData0.SL, TData0.FX, SL_vec', FX0_fz_nom_vec', ...
  '$\kappa [-]$', '$F_{x0}$ [N]', data_label, 'fig_fit_nominal_conditions_FX', ...
  'Fitting in nominal conditions', line_width, font_size_title)


%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
% extract data with variable load
[TDataDFz, ~] = intersect_table_data(SA_0, GAMMA_0);

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
%FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TDataDFz.SL));
ones_vec  = ones(size(TDataDFz.SL));

% Guess values for parameters to be optimised
%    [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
P0 = [  0,   0,   0,  0,   0,   0,   0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised 
%    [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
lb = []; %[-1 -1 -1 -1 -1 -1 -1];
ub = []; % [1 1 1 1 1 1 1];

KAPPA_vec = TDataDFz.SL;
FX_vec    = TDataDFz.FX;
FZ_vec    = TDataDFz.FZ;

% check guess
SL_vec = -0.3:0.001:0.3;

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varFz(P,FX_vec, KAPPA_vec, 0, FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);
res_Fx0_varFz = fval
P_dfz                        
RMS_Fx0_varFz = sqrt(res_Fx0_varFz*sum(FX_vec.^2)/length(FX_vec));
R2_Fx0_varFz = 1 - res_Fx0_varFz*sum(FX_vec.^2)/sum((FX_vec - mean(FX_vec)).^2);

% Change tyre data with new optimal values    
tyre_coeffs.pDx2 = P_dfz(1); % 1
tyre_coeffs.pEx2 = P_dfz(2);  
tyre_coeffs.pEx3 = P_dfz(3);
tyre_coeffs.pHx2 = P_dfz(4);
tyre_coeffs.pKx2 = P_dfz(5); 
tyre_coeffs.pKx3 = P_dfz(6);
tyre_coeffs.pVx2 = P_dfz(7);
%%
tmp_zeros = zeros(size(SL_vec));
tmp_ones = ones(size(SL_vec));

FX0_fz_var_vec1 = MF96_FX0_vec(SL_vec,tmp_zeros, tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec2 = MF96_FX0_vec(SL_vec,tmp_zeros, tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec3 = MF96_FX0_vec(SL_vec,tmp_zeros, tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec4 = MF96_FX0_vec(SL_vec,tmp_zeros, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

x_fit = repmat(SL_vec', 1, 4);
y_fit = [FX0_fz_var_vec1', FX0_fz_var_vec2', FX0_fz_var_vec3', FX0_fz_var_vec4'];
data_label = ["220 [N]", "700 [N]", "900 [N]", "1120 [N]"];
%figure 6
plot_fitted_data(TDataDFz.SL, TDataDFz.FX, x_fit, y_fit, '$\kappa$ [-]',...
  '$F_{x0}$ [N]', data_label,'fig_fit_variable_load_FX', ...
  'Fitting with variable load, $\gamma$ =0 [deg]', line_width, font_size_title)

%%
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);

Calfa_vec1 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

% figure 7
plot_stiffness(SL_vec,FZ_220, FZ_700, FZ_900, FZ_1120, Calfa_vec1_0, ...
  Calfa_vec2_0, Calfa_vec3_0, Calfa_vec4_0, Calfa_vec1, Calfa_vec2, ...
  Calfa_vec3, Calfa_vec4,'Longitudinal Stiffness','longitudinal_stiffness_FX',...
  font_size_title ) 

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER
%--------------------------------------------------------------------------
% extract data with variable load
[TDataGamma, ~] = intersect_table_data(SA_0, FZ_220);

% Fit the coeffs { pDx3}

% Guess values for parameters to be optimised
P0 = [0]; 

% NOTE: many local minima => limits on parameters are fundamentals
lb = []; %[0];
ub = []; %[20];

zeros_vec = zeros(size(TDataGamma.SL));
ones_vec  = ones(size(TDataGamma.SL));

KAPPA_vec = TDataGamma.SL;
GAMMA_vec = TDataGamma.IA; 
FX_vec    = TDataGamma.FX;
FZ_vec    = TDataGamma.FZ;

% NON CALNCELLARE
fig_camber_FX = figure('Color', 'w');
plot(KAPPA_vec,FX_vec, '.', 'DisplayName', '$F_{z}=220 [N]$');
grid on
xlabel('$\gamma$ [rad]')
ylabel('$F_X$ [N]')
title('Longitudianl force as function of camber angle')
legend('Location', 'best')
export_fig(fig_camber_FX, 'images\fig_camber_FX.png')

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varGamma(P,FX_vec, KAPPA_vec, GAMMA_vec, tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);
res_Fx0_varGamma = fval
P_varGamma
RMS_Fx0_varGamma = sqrt(res_Fx0_varGamma*sum(FX_vec.^2)/length(FX_vec));
R2_Fx0_varGamma = 1 - res_Fx0_varGamma*sum(FX_vec.^2)/sum((FX_vec - mean(FX_vec)).^2);

% Change tyre data with new optimal values
tyre_coeffs.pDx3 = P_varGamma(1); % 1

%%
FX0_varGamma_vec = MF96_FX0_vec(KAPPA_vec, zeros_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

TDataGamma0 = intersect_table_data(GAMMA_0, FZ_220);
TDataGamma2 = intersect_table_data(GAMMA_2, FZ_220);
TDataGamma4 = intersect_table_data(GAMMA_4, FZ_220);
FX0_Gamma0 = MF96_FX0_vec(TDataGamma0.SL, zeros(size(TDataGamma0.SL)), TDataGamma0.IA, tyre_coeffs.FZ0*ones(size(TDataGamma0.SL)), tyre_coeffs);
FX0_Gamma2 = MF96_FX0_vec(TDataGamma2.SL, zeros(size(TDataGamma2.SL)), TDataGamma2.IA, tyre_coeffs.FZ0*ones(size(TDataGamma2.SL)), tyre_coeffs);
FX0_Gamma4 = MF96_FX0_vec(TDataGamma4.SL,zeros(size(TDataGamma4.SL)), TDataGamma4.IA, tyre_coeffs.FZ0*ones(size(TDataGamma4.SL)), tyre_coeffs);
plot_y = cell(3,1);
plot_x = cell(3,1);
plot_fit = cell(3,1);
plot_x = {TDataGamma0.SL, TDataGamma2.SL, TDataGamma4.SL};
plot_y = {TDataGamma0.FX, TDataGamma2.FX, TDataGamma4.FX};
plot_fit = {FX0_Gamma0, FX0_Gamma2, FX0_Gamma4};
plot_label = ["$\gamma = 0 [deg]$", "$\gamma = 2 [deg]$", "$\gamma = 4 [deg]$"];

plot_fitted_data_struct(plot_x, plot_y, plot_x, plot_fit, ...
  '$\kappa$ [-]', '$F_{x0}$ [N]', plot_label, 'fig_fit_variable_camber_FX', ...
  'Fitting with variable camber', line_width, font_size_title, colors_vect);

data_label = ["Fitted data FZ=220 [N]"];
x_fit = repmat(KAPPA_vec, 1, 2);
y_fit = [FX0_varGamma_vec];
plot_fitted_data(TDataGamma.SL, TDataGamma.FX, x_fit, y_fit, ...
  '$\kappa$ [-]', '$F_{x0}$ [N]', data_label, ...
  'fig_fit_variable_camber_FX_1plot', 'Fitting with variable camber', ...
  line_width, font_size_title)

%%
fprintf('\n')
[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('mux     = %6.3f\n',Dx/tyre_coeffs.FZ0);
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);

% SSE is the sum of squared error,  SST is the sum of squared total
%% CHECK THE RESIDUALS
fprintf("Nominal conditions: %6.3f\n", res_Fx0);
fprintf("Variable load: %6.3f\n", res_Fx0_varFz);
fprintf("Variable camber: %6.3f\n", res_Fx0_varGamma);
fprintf('\n')
fprintf('R-squared = %6.3f\n',R2_Fx0);
fprintf('R-squared = %6.3f\n',R2_Fx0_varFz);
fprintf('R-squared = %6.3f\n',R2_Fx0_varGamma);
fprintf('\n')
fprintf('RMS = %6.3f\n', RMS_Fx0);
fprintf('RMS = %6.3f\n', RMS_Fx0_varFz);
fprintf('RMS = %6.3f\n', RMS_Fx0_varGamma);


%% Save tyre data structure to mat file
name = ["resFxzero", "resFxzerovarFz", "resFxzerovarGamma", "RtwoFxzero", ...
  "RtwoFxzerovarFz", "RtwoFxzerovarGamma", "RMSFxzero", "RMSFxzerovarFz", "RMSFxzerovarGamma"];
data = [res_Fx0, res_Fx0_varFz, res_Fx0_varGamma, R2_Fx0, R2_Fx0_varFz, R2_Fx0_varGamma, RMS_Fx0, RMS_Fx0_varFz, RMS_Fx0_varGamma];
write_latex_macro('results_to_lates_FX0.tex', name, data, 'w');
save(['tyre_' struct_name,'.mat'],'tyre_coeffs');



