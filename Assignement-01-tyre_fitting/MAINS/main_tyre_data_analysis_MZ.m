%% FITTING TYRE DATA-SELF ALIGNING MOMENT WITH MAGIC FORMULA 96

%% Initialisation
clc
clearvars 
close all 

addpath('utilities\')
addpath('tyre_lib\MZ')


warning('off', 'MATLAB:Figure:SetPosition')


% Choice of the dataset
data_set_path = 'dataset/';
data_set = 'Hoosier_B1464run23.mat';
struct_name = 'Hoosier_B1464';  
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
cut_start = 27760;%9000;
cut_end = 54500;%34500;%length(FY);
smpl_range = cut_start:cut_end;

% figure 1
plots_raw_data(cut_start,cut_end,FZ,IA,SA,SL,P,TSTC,TSTI,TSTO, ...
  font_size_title,'raw_data_MZ');

%% Sort data
sorting_data;

% figure 2
plot_sorted_data(tyre_data, idx, vec_samples, GAMMA_0, GAMMA_1, ...
  GAMMA_2, GAMMA_3, GAMMA_4, GAMMA_5, FZ_220, FZ_440, FZ_700, FZ_900, FZ_1120, ...
  FZ_1550, load_type, SA_0, SA_3neg, SA_6neg, font_size_title, ...
  'Sorted data self aligning','sorted_data_MZ0');

%% Intersect tables to obtain specific sub-datasets and plot them
[TData0, ~] = intersect_table_data(GAMMA_0, FZ_220 );

% figure 3
 plot_selected_data(TData0, font_size_title);

%%-------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
if exist(['tyre_' struct_name,'.mat'], 'file')
  load(['tyre_' struct_name,'.mat']);
  tyre_coeffs.qHz1 = 0; 
  tyre_coeffs.qBz1 = 0;
  tyre_coeffs.qCz1 = 0;
  tyre_coeffs.qDz1 = 0;
  tyre_coeffs.qEz1 = 0;
  tyre_coeffs.qEz4 = 0;
  tyre_coeffs.qBz9 = 0;
  tyre_coeffs.Bz10  = 0;
  tyre_coeffs.qDz6 = 0;
  tyre_coeffs.qHz2 = 0;
  tyre_coeffs.qBz2 = 0;
  tyre_coeffs.qBz3 = 0;
  tyre_coeffs.qDz2 = 0;
  tyre_coeffs.qEz2 = 0;
  tyre_coeffs.qEz3 = 0;
  tyre_coeffs.qDz7 = 0;
  tyre_coeffs.qHz3 = 0;
  tyre_coeffs.qBz4 = 0;
  tyre_coeffs.qBz5 = 0;
  tyre_coeffs.qDz3 = 0;
  tyre_coeffs.qDz4 = 0;
  tyre_coeffs.qEz5 = 0;
  tyre_coeffs.qDz8 = 0;
  tyre_coeffs.qHz4 = 0;
  tyre_coeffs.qDz9 = 0;
   
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

%% ------------------------------------------------------------------------
% FITTING WITH NOMINAL LOAD Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10
%--------------------------------------------------------------------------
% Guess values for parameters to be optimised
%    [qHz1, qBz1, qCz1, qDz1, qEz1, qEz4, qBz9, qDz6, qBz10]
P0 = -0.75*ones(1, 9);

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [qHz1, qBz1, qCz1, qDz1, qEz1, qEz4, qBz9, qDz6, qBz10]
lb = []; % lower bound
ub = []; % upper bound

ALPHA_vec = TData0.SA;  % lateral slip angle
MZ_vec    = TData0.MZ;  % aligning moment
FY_vec    = TData0.FY;  % laterla force

% check guess
SA_vec = -0.3:0.001:0.3; % lateral slip to be used in the plot for check 
% the interpolation outside the input data range

% resid_pure_Mz returns the residual, so minimize the residual varying Y. 
% It is an unconstrained minimization problem 
[P_fz_nom, res_Mz0, exitflag] = fmincon(@(P)resid_pure_Mz(P, MZ_vec, ...
  ALPHA_vec, 0, FZ0, FY_vec, tyre_coeffs),P0,[],[],[],[],lb,ub);

P_fz_nom
res_Mz0
RMS_Mz0 = sqrt(res_Mz0*sum(MZ_vec.^2)/length(MZ_vec));
R2_Mz0 = 1 - res_Mz0*sum(MZ_vec.^2)/sum((MZ_vec - mean(MZ_vec)).^2);

% Update tyre data with new optimal values                   
tyre_coeffs.qHz1  = P_fz_nom(1); 
tyre_coeffs.qBz1  = P_fz_nom(2);
tyre_coeffs.qCz1  = P_fz_nom(3);
tyre_coeffs.qDz1  = P_fz_nom(4);
tyre_coeffs.qEz1  = P_fz_nom(5);
tyre_coeffs.qEz4  = P_fz_nom(6);
tyre_coeffs.qBz9  = P_fz_nom(7);
tyre_coeffs.qDz6  = P_fz_nom(8);
tyre_coeffs.qBz10 = P_fz_nom(9);

MZ0_fz_nom_vec = MF96_MZ0_vec(SA_vec, zeros(size(SA_vec)), FZ0.*ones(size(SA_vec)), tyre_coeffs);
        
data_label = "Fitted data"; 
% figure 5          
plot_fitted_data(TData0.SA, TData0.MZ, SA_vec', MZ0_fz_nom_vec', ...
  '$\alpha [-]$', '$M_{z0}$ [Nm]', data_label, ...
  'fig_fit_nominal_conditions_MZ0', 'Fitting in nominal conditions', ...
  line_width, font_size_title)


%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
% extract data with variable load
[TDataDFz, ~] = GAMMA_0; % intersect_table_data(SA_0, GAMMA_0);

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
%FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [ qHz2, qBz2, qBz3, qDz2, qEz2, qEz3, qDz7]
P0 = 0*ones(7,1); % [0, 0, 0, 0, 0, 0, 0];

lb = []; %[0, 0, 0, 0, 0, 0, 0];
ub = []; %[10, 10, 10, 10, 10, 10, 10];

ALPHA_vec = TDataDFz.SA;
MZ_vec    = TDataDFz.MZ;
FZ_vec    = TDataDFz.FZ;

% check guess
SA_vec = -0.3:0.001:0.3;
% FX0_dfz_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)), zeros(size(SL_vec)), ...
%                            TDataDFz.FZ,tyre_coeffs);

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_dfz,res_Mz0_varFz,exitflag] = fmincon(@(P)resid_pure_Mz_varFz(P, MZ_vec, ...
  ALPHA_vec, 0, FZ_vec, tyre_coeffs),P0,[],[],[],[],lb,ub);
P_dfz
res_Mz0_varFz
RMS_Mz0_varFz = sqrt(res_Mz0_varFz*sum(MZ_vec.^2)/length(MZ_vec));
R2_Mz0_varFz = 1 - res_Mz0_varFz*sum(MZ_vec.^2)/sum((MZ_vec - mean(MZ_vec)).^2);

% Change tyre data with new optimal values          
tyre_coeffs.qHz2 = P_dfz(1); 
tyre_coeffs.qBz2 = P_dfz(2); 
tyre_coeffs.qBz3 = P_dfz(3); 
tyre_coeffs.qDz2 = P_dfz(4); 
tyre_coeffs.qEz2 = P_dfz(5); 
tyre_coeffs.qEz3 = P_dfz(6); 
tyre_coeffs.qDz7 = P_dfz(7); 

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));

%%
MZ0_fz_var_vec1 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec2 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec3 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec4 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec5 = MF96_MZ0_vec(SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones, tyre_coeffs);

x_fit = repmat(SA_vec', 1, 5);
y_fit = [MZ0_fz_var_vec1', MZ0_fz_var_vec2', MZ0_fz_var_vec3', MZ0_fz_var_vec4', MZ0_fz_var_vec5'];
data_label = string(5);
data_label = ["$F_{z}$ = 220 [N]", "$F_{z}$ = 440 [N]", ...
  "$F_{z}$ = 700 [N]", "$F_{z}$ = 900 [N]", "$F_{z}$ = 1120 [N]"];
%figure 6
plot_fitted_data(TDataDFz.SA, TDataDFz.MZ, x_fit, y_fit, '$\alpha$ [-]',...
  '$M_{z0}$ [Nm]', data_label,'fig_fit_variable_load_MZ0', ...
  'Fitting with variable load', line_width, font_size_title)

T220  = intersect_table_data(FZ_220, GAMMA_0);
T440  = intersect_table_data(FZ_440, GAMMA_0);
T700  = intersect_table_data(FZ_700, GAMMA_0);
T900  = intersect_table_data(FZ_900, GAMMA_0);
T1120 = intersect_table_data(FZ_1120, GAMMA_0);

SA_cell = cell(1, 5);
SA_cell = {T220.SA, T440.SA, T700.SA, T900.SA, T1120.SA};
MZ_cell = cell(1, 5);
MZ_cell = {T220.MZ, T440.MZ, T700.MZ, T900.MZ, T1120.MZ};
x_fit_cell = cell(1, 5);
y_fit_cell = cell(1, 5);
for i=1:5
  x_fit_cell{i} = x_fit(:, i);
  y_fit_cell{i} = y_fit(:, i);
end
plot_fitted_data_struct(SA_cell, MZ_cell, x_fit_cell, y_fit_cell, '$\alpha$ [-]',...
  '$M_{z0}$ [Nm]', data_label,'fig_fit_variable_load_MZ0_multicolot', ...
  'Fitting with variable load', line_width, font_size_title, colors_vect);

%% Stiffness
[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_220.FZ), tyre_coeffs);
[~, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, By, Cy, Dy, Ey, SVy, mean(FZ_220.FZ), tyre_coeffs);

[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_440.FZ), tyre_coeffs);
[~, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_440.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, By, Cy, Dy, Ey, SVy, mean(FZ_440.FZ), tyre_coeffs);

[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_700.FZ), tyre_coeffs);
[~, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, By, Cy, Dy, Ey, SVy, mean(FZ_700.FZ), tyre_coeffs);

[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_900.FZ), tyre_coeffs);
[~, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, By, Cy, Dy, Ey, SVy, mean(FZ_900.FZ), tyre_coeffs);

[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_1120.FZ), tyre_coeffs);
[~, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec5_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, By, Cy, Dy, Ey, SVy, mean(FZ_1120.FZ), tyre_coeffs);





% [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_440.FZ), tyre_coeffs);
% Calfa_vec2_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, mean(FZ_440.FZ), tyre_coeffs);
% [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_700.FZ), tyre_coeffs);
% Calfa_vec3_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, mean(FZ_700.FZ), tyre_coeffs);
% [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_900.FZ), tyre_coeffs);
% Calfa_vec4_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, mean(FZ_900.FZ), tyre_coeffs);
% [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_1120.FZ), tyre_coeffs);
% Calfa_vec5_0 = magic_formula_stiffness_MZ(0, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, mean(FZ_1120.FZ), tyre_coeffs);

Calfa_vec1 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_220.FZ)  * tmp_ones, tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_440.FZ)  * tmp_ones, tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_700.FZ)  * tmp_ones, tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_900.FZ)  * tmp_ones, tyre_coeffs);
Calfa_vec5 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_1120.FZ) * tmp_ones, tyre_coeffs);

%% figure 7
plot_stiffness_MZ; 
% plot_stiffness_FY(SA_vec,FZ_220, FZ_700, FZ_900, FZ_1120, FZ_1550, ...
%   Calfa_vec1_0,Calfa_vec2_0, Calfa_vec3_0, Calfa_vec4_0, Calfa_vec5_0, ...
%   Calfa_vec1, Calfa_vec2, Calfa_vec3, Calfa_vec4, Calfa_vec5, ...
%   'Cornering stiffness', 'cornering_stiffness_MZ', font_size_title )
%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER
%--------------------------------------------------------------------------
% extract data with variable load
[TDataGamma, ~] = FZ_220; % intersect_table_data(SA_0, FZ_220);

% Guess values for parameters to be optimised
%    [qHz3, qBz4, qBz5, qDz3, qDz4, qEz5, qDz8, qHz4, qDz9]
P0 = 1*ones(1, 9); %[1, 1, 1, 1, 1, 1, 1, 1, 1]; 

% NOTE: many local minima => limits on parameters are fundamentals
lb = [];
ub = [];

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

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma, res_Mz0_varGamma, exitflag] = fmincon(@(P)resid_pure_Mz_varGamma(P, MZ_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0, FY_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);
P_varGamma
res_Mz0_varGamma
RMS_Mz0_varGamma = sqrt(res_Mz0_varGamma*sum(MZ_vec.^2)/length(MZ_vec));
R2_Mz0_varGamma = 1 - res_Mz0_varGamma*sum(MZ_vec.^2)/sum((MZ_vec - mean(MZ_vec)).^2);

% Change tyre data with new optimal values
tyre_coeffs.qHz3 = P_varGamma(1); 
tyre_coeffs.qBz4 = P_varGamma(2); 
tyre_coeffs.qBz5 = P_varGamma(3); 
tyre_coeffs.qDz3 = P_varGamma(4); 
tyre_coeffs.qDz4 = P_varGamma(5); 
tyre_coeffs.qEz5 = P_varGamma(6); 
tyre_coeffs.qDz8 = P_varGamma(7); 
tyre_coeffs.qHz4 = P_varGamma(8); 
tyre_coeffs.qDz9 = P_varGamma(9);
%%
alpha_vec = -0.25:0.001:0.25;
ones_vec = ones(length(alpha_vec), 1);
MZ0_Gamma0 = MF96_MZ0_vec(alpha_vec, 0*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
MZ0_Gamma1 = MF96_MZ0_vec(alpha_vec, 1*pi/180*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
MZ0_Gamma2 = MF96_MZ0_vec(alpha_vec, 2*pi/180*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
MZ0_Gamma3 = MF96_MZ0_vec(alpha_vec, 3*pi/180*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
MZ0_Gamma4 = MF96_MZ0_vec(alpha_vec, 4*pi/180*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);

fig_fit_variable_camber_MZ = figure('Color', 'w');
hold on
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
grid on
export_fig(fig_fit_variable_camber_MZ, 'images\fig_fit_variable_camber_MZ0.png');
%%
x_fit_cell = cell(1, 5);
y_fit_cell = cell(1, 5);
y_raw_cell = cell(1, 5);
x_raw_cell = cell(1, 5);
x_raw_cell = {TDataGamma0.SA, TDataGamma1.SA, TDataGamma2.SA, TDataGamma3.SA, TDataGamma4.SA};
x_fit_cell = {alpha_vec, alpha_vec, alpha_vec, alpha_vec, alpha_vec};
y_fit_cell = {MZ0_Gamma0, MZ0_Gamma1, MZ0_Gamma2, MZ0_Gamma3, MZ0_Gamma4};
y_raw_cell = {TDataGamma0.MZ, TDataGamma1.MZ, TDataGamma2.MZ, TDataGamma3.MZ, TDataGamma4.MZ};
data_label = ["$\gamma = 0 [deg]$", "$\gamma = 1 [deg]$", "$\gamma = 2 [deg]$", "$\gamma = 3 [deg]$", "$\gamma = 4 [deg]$"];
plot_fitted_data_struct(x_raw_cell, y_raw_cell, x_fit_cell, y_fit_cell, ...
  '$\alpha []$', '$M_{z0}$ [Nm]', data_label,  'fig_fit_variable_camber_MZ0_2', ...
  'Self alignign moment for variable camber', line_width, font_size_title, colors_vect);


% NON CANCELLARE
% fig_camber_MZ = figure('Color', 'w');
% plot(ALPHA_vec, MZ_vec, '.', 'DisplayName', '$F_{z}=220 [N]$');
% grid on
% xlabel('$\gamma$ [rad]')
% ylabel('$M_Z$ [N]')
% title('Self aligning moment as function of camber angle')
% legend('Location', 'best')
% export_fig(fig_camber_MZ, 'images\fig_camber_MZ.png')

%%
ones_vec = ones(length(ALPHA_vec), 1);
MZ0_varGamma_vec = MF96_MZ0_vec(ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

data_label = ["Fitted data FZ=220 [N]"];
x_fit = repmat(ALPHA_vec, 1, 2);
y_fit = [MZ0_varGamma_vec];
plot_fitted_data(TDataGamma.SA, TDataGamma.MZ, x_fit, y_fit, '$\alpha$ [-]', ...
  '$M_{z0}$ [N]', data_label, 'fig_fit_variable_camber_MZ0', ...
  'Fitting with variable camber', line_width, font_size_title)


fprintf("Nominal conditions: %6.3f\n", res_Mz0);
fprintf("Variable load: %6.3f\n", res_Mz0_varFz);
fprintf("Variable camber: %6.3f\n", res_Mz0_varGamma);
fprintf('\n')
fprintf('R-squared = %6.3f\n',R2_Mz0);
fprintf('R-squared = %6.3f\n',R2_Mz0_varFz);
fprintf('R-squared = %6.3f\n',R2_Mz0_varGamma);
fprintf('\n')
fprintf('RMS = %6.3f\n', RMS_Mz0);
fprintf('RMS = %6.3f\n', RMS_Mz0_varFz);
fprintf('RMS = %6.3f\n', RMS_Mz0_varGamma); 


%% Save tyre data structure to mat file
name = ["resMzzero", "resMzzerovarFz", "resMzzerovarGamma", "RtwoMzzero", ...
  "RtwoMzzerovarFz", "RtwoMzzerovarGamma", "RMSMzzero", "RMSMzzerovarFz", "RMSMzzerovarGamma"];
data = [res_Mz0, res_Mz0_varFz, res_Mz0_varGamma, R2_Mz0, R2_Mz0_varFz, R2_Mz0_varGamma, RMS_Mz0, RMS_Mz0_varFz, RMS_Mz0_varGamma];
write_latex_macro('results_to_lates_MZ0.tex', name, data, 'w');


%% Save tyre data structure to mat file
%
save(['tyre_' struct_name,'.mat'],'tyre_coeffs');