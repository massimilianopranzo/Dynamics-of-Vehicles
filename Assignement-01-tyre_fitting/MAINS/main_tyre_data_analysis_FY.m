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

% addpath('..\utilities\')
% addpath('..\dataset\')
addpath('\MAINS\..') 
addpath('tyre_lib\FY')

% Suppress warning export_fig for docked figures
warning('off', 'MATLAB:Figure:SetPosition')

% Choice of the dataset
data_set_path = 'dataset/';
data_set = 'Hoosier_B1464run23.mat';
struct_name = 'Hoosier_B1464';  
load_type = 'lateral'; 

initialization


%% Load raw data  
load ([data_set_path, data_set]); % pure lateral

% select dataset portion
cut_start = 27760;%9000;
cut_end = 54500;%34500;%length(FY);
smpl_range = cut_start:cut_end;

% figure 1
plots_raw_data(cut_start,cut_end,FZ,IA,SA,SL,P,TSTC,TSTI,TSTO, ...
  font_size_title,'raw_data_FY');

%% Sort data
sorting_data;

% figure 2
plot_sorted_data(tyre_data, idx, vec_samples, GAMMA_0, GAMMA_1, ...
  GAMMA_2, GAMMA_3, GAMMA_4, GAMMA_5, FZ_220, FZ_440, FZ_700, FZ_900, FZ_1120, ...
  FZ_1550, load_type, SA_0, SA_3neg, SA_6neg, font_size_title, ...
  'Sorted data lateral','sorted_data_FY');

%% -------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
[TData0, ~] = intersect_table_data(GAMMA_0, FZ_220);

% figure 3
plot_selected_data(TData0, font_size_title);
if exist(['tyre_', struct_name,'.mat'], 'file')
  load(['tyre_', struct_name,'.mat']);
  tyre_coeffs.pHy1 = 0;
  tyre_coeffs.pCy1 = 0;
  tyre_coeffs.pDy1 = 0;
  tyre_coeffs.pKy1 = 0;
  tyre_coeffs.pKy2 = 0;
  tyre_coeffs.pEy1 = 0;
  tyre_coeffs.pVy1 = 0;
  tyre_coeffs.pHy2 = 0;
  tyre_coeffs.pDy2 = 0;
  tyre_coeffs.pEy2 = 0;
  tyre_coeffs.pVy2 = 0;
  tyre_coeffs.pHy3 = 0;
  tyre_coeffs.pDy3 = 0;
  tyre_coeffs.pKy3 = 0;
  tyre_coeffs.pEy3 = 0;
  tyre_coeffs.pEy4 = 0;
  tyre_coeffs.pVy3 = 0;
  tyre_coeffs.pVy4 = 0;
else
  tyre_coeffs = initialise_tyre_data(R0, Fz0);
end

FZ0 = mean(TData0.FZ);
zeros_vec = zeros(size(TData0.SA));
ones_vec  = ones(size(TData0.SA));

%% ------------------------------------------------------------------------
% FITTING WITH NOMINAL LOAD Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10
%--------------------------------------------------------------------------
% Guess values for parameters to be optimised
%    [pCy1 pDy1 pEy1  pHy1  pKy1  pKy2, pVy1]
P1 = [  0.7,   3.1,   -0.4,     0,   -110,     -3.2,  0]; % parametri stefano
% P0 = [1 1 1 1 1 1 1]; 
name1 = ["pCy1", "pDy1", "pEy1", "pHy1", "pKy1", "pKy2", "pVy1"];

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [pCy1 pDy1 pEy1  pHy1  pKy1  pKy2, pVy1]
lb1 = [1   -inf -inf -inf -inf -inf -inf]; % lower bound -inf*ones(1, 7);
ub1 = [inf  inf  1    inf  inf  inf  inf];%[2 inf inf inf inf inf inf]; % upper bound % inf*ones(1, 7);

ALPHA_vec = TData0.SA;  % side slip angle
FY_vec    = TData0.FY;  % lateral force

% check guess
SA_vec = -0.28:0.001:0.28; % lateral slip to be used in the plot for check the interpolation outside the input data range

% resid_pure_Fy returns the residual, so minimize the residual varying Y. 
% It is an unconstrained minimization problem 
[P_fz_nom, res_Fy0, exitflag] = fmincon(@(P)resid_pure_Fy(P, FY_vec, ALPHA_vec, 0, FZ0, tyre_coeffs),...
                               P1,[],[],[],[],lb1,ub1);
P_fz_nom
res_Fy0
RMS_Fy0 = sqrt(res_Fy0*sum(FY_vec.^2)/length(FY_vec));
R2_Fy0 = 1 - res_Fy0*sum(FY_vec.^2)/sum((FY_vec - mean(FY_vec)).^2);

% Update tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDy1 = P_fz_nom(2) ;  
tyre_coeffs.pEy1 = P_fz_nom(3) ;
tyre_coeffs.pHy1 = P_fz_nom(4) ; 
tyre_coeffs.pKy1 = P_fz_nom(5) ;
tyre_coeffs.pKy2 = P_fz_nom(6) ;
tyre_coeffs.pVy1 = P_fz_nom(7) ;

FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs);
       
data_label = "Fitted data"; 
% figure 5          
plot_fitted_data(TData0.SA, TData0.FY, SA_vec', FY0_fz_nom_vec', ...
  '$\alpha [-]$', '$F_{y0}$ [N]', data_label, ...
  'fig_fit_nominal_conditions_FY', 'Fitting in nominal conditions', line_width, ...
  font_size_title)


%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
% extract data with variable load
[TDataDFz, ~] = GAMMA_0;

zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));

% Guess values for parameters to be optimised
%    [   pDy2     pEy2   pHy2     pVy2    ]

P2 = [0 0 0 0];
name2 = ["pDy2", "pEy2", "pHy2", "pVy2"];

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [   pDy2     pEy2   pHy2     pVy2  ]
% lb = [-5 -20 -5 -5] ;
% ub = [15 30 5 10];
lb2 = [];
ub2 = [];


ALPHA_vec = TDataDFz.SA;
FY_vec    = TDataDFz.FY;
FZ_vec    = TDataDFz.FZ;

% check guess
SA_vec = -0.3:0.001:0.3;

% LSM_pure_Fy returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_dfz, res_Fy0_varFz, exitflag] = fmincon(@(P)resid_pure_Fy_varFz(P, ...
  FY_vec, ALPHA_vec, 0, FZ_vec, tyre_coeffs),P2,[],[],[],[],lb2,ub2);
P_dfz
res_Fy0_varFz
RMS_Fy0_varFz = sqrt(res_Fy0_varFz*sum(FY_vec.^2)/length(FY_vec));
R2_Fy0_varFz = 1 - res_Fy0_varFz*sum(FY_vec.^2)/sum((FY_vec - mean(FY_vec)).^2);

% Change tyre data with new optimal values                             
tyre_coeffs.pDy2 = P_dfz(1); % 1
tyre_coeffs.pEy2 = P_dfz(2);  
tyre_coeffs.pHy2 = P_dfz(3);
tyre_coeffs.pVy2 = P_dfz(4);

% figure
plot_variable_loads

%% figure 7
% plot_stiffness_FY; 
plot_stiffness_FY(SA_vec,FZ_220, FZ_440, FZ_700, FZ_900, FZ_1120, ...
  Calfa_vec1_0,Calfa_vec2_0, Calfa_vec3_0, Calfa_vec4_0, Calfa_vec5_0, ...
  Calfa_vec1, Calfa_vec2, Calfa_vec3, Calfa_vec4, Calfa_vec5, ...
  'Cornering stiffness', 'cornering_stiffness_FY', font_size_title);

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER
%--------------------------------------------------------------------------
% extract data with variable load
[TDataGamma, ~] = FZ_220; % intersect_table_data(SA_0, FZ_220);

% Guess values for parameters to be optimised
%    [pHy3, pDy3, pKy3, pEy3, pEy4, pVy3, pVy4]
P3 = [ 0     5     2     0.9     -5     0    -4] ; 
name3 = ["pHy3", "pDy3", "pKy3", "pEy3", "pEy4", "pVy3", "pVy4"];
% NOTE: many local minima => limits on parameters are fundamentals
lb3 = [ ];
ub3 = [ ];

zeros_vec = zeros(size(TDataGamma.SA));
ones_vec  = ones(size(TDataGamma.SA));

ALPHA_vec = TDataGamma.SA;
GAMMA_vec = TDataGamma.IA; 
FY_vec    = TDataGamma.FY;
FZ_vec    = TDataGamma.FZ;

% NOT CANCELLARE
fig_camber_FY = figure('Color', 'w');
plot(ALPHA_vec,FY_vec, '.', 'DisplayName', '$F_{z}=220 [N]$');
grid on
xlabel('$\gamma$ [rad]')
ylabel('$F_{y0}$ [N]')
title('Lateral force as function of camber angle')
legend('Location', 'best')
export_fig(fig_camber_FY, 'images\fig_camber_FY.png')

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma, res_Fy0_varGamma, exitflag] = fmincon(@(P) ...
  resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec, GAMMA_vec, ...
  tyre_coeffs.FZ0*ones_vec , tyre_coeffs),P3,[],[],[],[],lb3, ...
  ub3); %tyre_coeffs.FZ0,
P_varGamma
res_Fy0_varGamma
RMS_Fy0_varGamma = sqrt(res_Fy0_varGamma*sum(FY_vec.^2)/length(FY_vec));
R2_Fy0_varGamma = 1 - res_Fy0_varGamma*sum(FY_vec.^2)/sum((FY_vec - mean(FY_vec)).^2);

%%
% Change tyre data with new optimal values
tyre_coeffs.pHy3 = P_varGamma(1); 
tyre_coeffs.pDy3 = P_varGamma(2); 
tyre_coeffs.pKy3 = P_varGamma(3); 
tyre_coeffs.pEy3 = P_varGamma(4); 
tyre_coeffs.pEy4 = P_varGamma(5); 
tyre_coeffs.pVy3 = P_varGamma(6); 
tyre_coeffs.pVy4 = P_varGamma(7); 

FY0_varGamma_vec = MF96_FY0_vec(zeros_vec, ALPHA_vec, GAMMA_vec, ...
  tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% data_label = ["Fitted data FZ=220 [N]"];
% x_fit = repmat(ALPHA_vec, 1, 2);
% y_fit = [FY0_varGamma_vec];
% plot_fitted_data(TDataGamma.SA, TDataGamma.FY, x_fit, y_fit, ...
%   '$\alpha$ [-]', '$F_{Y}$ [N]', data_label, 'fig_fit_variable_camber_FY', ...
%   'Fitting with variable camber', line_width, font_size_title)

%%
alpha_vec = -0.25:0.001:0.25;
ones_vec = ones(length(alpha_vec), 1 );
TDataGamma0 = intersect_table_data(GAMMA_0, FZ_220);
TDataGamma1 = intersect_table_data(GAMMA_1, FZ_220);
TDataGamma2 = intersect_table_data(GAMMA_2, FZ_220);
TDataGamma3 = intersect_table_data(GAMMA_3, FZ_220);
TDataGamma4 = intersect_table_data(GAMMA_4, FZ_220);
FY0_Gamma0 = MF96_FY0_vec(zeros_vec, alpha_vec, 0*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
FY0_Gamma1 = MF96_FY0_vec(zeros_vec, alpha_vec, -1*pi/180*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
FY0_Gamma2 = MF96_FY0_vec(zeros_vec, alpha_vec, -2*pi/180*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
FY0_Gamma3 = MF96_FY0_vec(zeros_vec, alpha_vec, -3*pi/180*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
FY0_Gamma4 = MF96_FY0_vec(zeros_vec, alpha_vec, -4*pi/180*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
%%
plot_y = cell(5,1);
plot_x = cell(5,1);
plot_fit = cell(5,1);
x_fit = cell(5, 1);
plot_x = {TDataGamma0.SA, TDataGamma1.SA, TDataGamma2.SA, TDataGamma3.SA, TDataGamma4.SA};
plot_y = {TDataGamma0.FY, TDataGamma1.FY, TDataGamma2.FY, TDataGamma3.FY, TDataGamma4.FY};
plot_fit = {FY0_Gamma0, FY0_Gamma1, FY0_Gamma2, FY0_Gamma3, FY0_Gamma4};
x_fit = {alpha_vec, alpha_vec, alpha_vec, alpha_vec, alpha_vec};
plot_label = ["$\gamma = 0 [deg]$", "$\gamma = 1 [deg]$", "$\gamma = 2 [deg]$", "$\gamma = 3 [deg]$", "$\gamma = 4 [deg]$"];

plot_fitted_data_struct(plot_x, plot_y, x_fit, plot_fit, ...
  '$\alpha$ [-]', '$F_{y0}$ [N]', plot_label, 'fig_fit_variable_camber_FY', ...
  'Fitting with variable camber - Nominal load', line_width, ...
  font_size_title, colors_vect);


%% -------------------
% PLOT THE RESIDUALS
% SSE is the sum of squared error,  SST is the sum of squared total
[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
fprintf('By      = %6.3f\n', By);
fprintf('Cy      = %6.3f\n', Cy);
fprintf('muy     = %6.3f\n', Dy/tyre_coeffs.FZ0);
fprintf('Ey      = %6.3f\n', Ey);
fprintf('SVy     = %6.3f\n', SVy);
fprintf('alpha_y = %6.3f\n', alpha__y);
fprintf('Ky      = %6.3f\n', By*Cy*Dy/tyre_coeffs.FZ0);


fprintf("Nominal conditions: %6.3f\n", res_Fy0);
fprintf("Variable load: %6.3f\n", res_Fy0_varFz);
fprintf("Variable camber: %6.3f\n", res_Fy0_varGamma);
fprintf('\n')
fprintf('R-squared = %6.3f\n',R2_Fy0);
fprintf('R-squared = %6.3f\n',R2_Fy0_varFz);
fprintf('R-squared = %6.3f\n',R2_Fy0_varGamma);
fprintf('\n')
fprintf('RMS = %6.3f\n', RMS_Fy0);
fprintf('RMS = %6.3f\n', RMS_Fy0_varFz);
fprintf('RMS = %6.3f\n', RMS_Fy0_varGamma); 

%%
P0 = {P1 P2 P3};
lb = {lb1 lb2 lb3};
ub = {ub1 ub2 ub3};
name = {name1 name2 name3};
table_P0("pure_lateral", name, P0, lb, ub)

%% Save tyre data structure to mat file
name = ["resFyzero", "resFyzerovarFz", "resFyzerovarGamma", "RtwoFyzero",...
  "RtwoFyzerovarFz", "RtwoFyzerovarGamma", "RMSFyzero", "RMSFyzerovarFz", "RMSFyzerovarGamma"];
data = [res_Fy0, res_Fy0_varFz, res_Fy0_varGamma, R2_Fy0, R2_Fy0_varFz, R2_Fy0_varGamma, RMS_Fy0, RMS_Fy0_varFz, RMS_Fy0_varGamma];
write_latex_macro('results_to_lates_FY0.tex', name, data, 'w');

save(['tyre_' struct_name,'.mat'],'tyre_coeffs');




