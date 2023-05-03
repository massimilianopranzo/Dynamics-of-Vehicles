%% Initialisation
clc
clearvars 
close all 

addpath('utilities\')
addpath('MAINS\..')
addpath('tyre_lib\')
addpath('tyre_lib\COMBINED\')

warning('off', 'MATLAB:Figure:SetPosition')

% Choice of the dataset
data_set_path = 'dataset/';
data_set = 'Hoosier_B1464run30.mat';
struct_name = 'Hoosier_B1464';  
load_type = 'longitudinal';

initialization

%% INITIALIZZARE

%% Load raw data  
load ([data_set_path, data_set]); % pure lateral

% select dataset portion
cut_start = 19028;
cut_end = 37643;
smpl_range = cut_start:cut_end;

% figure 1
plots_raw_data(cut_start,cut_end,FZ,IA,SA,SL,P,TSTC,TSTI,TSTO, ...
  font_size_title,'raw_data_longitudinal_FX');

%% Sort data
sorting_data;
% figure 2
plot_sorted_data(tyre_data, idx, vec_samples, GAMMA_0, GAMMA_1, ...
  GAMMA_2, GAMMA_3, GAMMA_4, GAMMA_5, FZ_220, FZ_440, FZ_700, FZ_900, FZ_1120, ...
  FZ_1550, load_type, SA_0, SA_3neg, SA_6neg, font_size_title, ...
  'Sorted data combined longitudinal','sorted_data_combined_FX');

if exist(['tyre_' struct_name,'.mat'], 'file')
  load(['tyre_' struct_name,'.mat']);
  % combined longitudinal
  tyre_coeffs.rHx1 = 0;
  tyre_coeffs.rBx1 = 0;
  tyre_coeffs.rCx1 = 0;
  tyre_coeffs.rBx2 = 0;

  % combined lateral
  tyre_coeffs.rVy1 = 0;
  tyre_coeffs.rVy4 = 0;
  tyre_coeffs.rVy5 = 0;
  tyre_coeffs.rVy6 = 0;
  tyre_coeffs.rHy1 = 0;
  tyre_coeffs.rBy1 = 0;
  tyre_coeffs.rBy2 = 0;
  tyre_coeffs.rBy3 = 0;
  tyre_coeffs.rCy1 = 0;
  tyre_coeffs.rVy2 = 0;
  tyre_coeffs.rVy3 = 0;
else
  error('Run the genertion of FX, FY and MZ parameters first!');
end


%% COMBINED LONGITUDINAL
% rHx1, rBx1, rCx1 rBx2
% -27.1725    5.5548   11.0463   -1.4499
P1 = [-27 5.5 11 -1.5];
lb1 = [-30 4 10 -2];
ub1 = [-20 7 15 0];
name1 = ["rHx1", "rBx1", "rCx1", "rBx2"];

FX_vec = tyre_data.FX;
KAPPA_vec = tyre_data.SL;
ALPHA_vec = tyre_data.SA;
GAMMA_vec = tyre_data.IA;
FZ_vec = tyre_data.FZ;

[P_min,res_Fx,exitflag] = fmincon(@(P)resid_long(P, FX_vec, KAPPA_vec, ALPHA_vec, GAMMA_vec, FZ_vec, tyre_coeffs),...
P1,[],[],[],[],lb1,ub1);

P_min
res_Fx
RMS_Fx = sqrt(res_Fx*sum(FX_vec.^2)/length(FX_vec));
R2_Fx = 1 - res_Fx*sum(FX_vec.^2)/sum((FX_vec - mean(FX_vec)).^2);

tyre_coeffs.rHx1 = P_min(1);
tyre_coeffs.rBx1 = P_min(2);
tyre_coeffs.rCx1 = P_min(3);
tyre_coeffs.rBx2 = P_min(4);      

FX_fit = cell(4, 1); % {FZ=220, FZ=700, FZ=900, FZ=1120}
FX_raw = cell(4, 1); % {FZ=220, FZ=700, FZ=900, FZ=1120}

Data_220 = intersect_table_data(FZ_220, GAMMA_0);
Data_700 = intersect_table_data(FZ_700, GAMMA_0);
Data_900 = intersect_table_data(FZ_900, GAMMA_0);
Data_1120 = intersect_table_data(FZ_1120, GAMMA_0);

[FX_fit{1}] = MF96_FXcomb(Data_220.SL, Data_220.SA, Data_220.IA, Data_220.FZ, tyre_coeffs);
[FX_fit{2}] = MF96_FXcomb(Data_700.SL, Data_700.SA, Data_700.IA, Data_700.FZ, tyre_coeffs);
[FX_fit{3}] = MF96_FXcomb(Data_900.SL, Data_900.SA, Data_900.IA, Data_900.FZ, tyre_coeffs);
[FX_fit{4}] = MF96_FXcomb(Data_1120.SL, Data_1120.SA, Data_1120.IA, Data_1120.FZ, tyre_coeffs);
FX_raw = {Data_220.FX, Data_700.FX, Data_900.FX, Data_1120.FX};
x_raw = {Data_220.SL, Data_700.SL, Data_900.SL, Data_1120.SL};


data_label = ["FZ=220 [N]", "FZ=700 [N]", "FZ=900 [N]", "FZ=1120 [N]"];
leg_angle = ["$\alpha$=0", "$\alpha=-3$", "$\alpha=-6$"];
%%
x_raw = cell(4, 1);
x_raw{1}{1}(:) = intersect_table_data(Data_220, SA_0).SL; % fz220
x_raw{1}{2}(:) = intersect_table_data(Data_220, SA_3neg).SL; % fz220
x_raw{1}{3}(:) = intersect_table_data(Data_220, SA_6neg).SL; % fz220
x_raw{2}{1}(:) = intersect_table_data(Data_700, SA_0).SL; % fz700
x_raw{2}{2}(:) = intersect_table_data(Data_700, SA_3neg).SL; % fz700
x_raw{2}{3}(:) = intersect_table_data(Data_700, SA_6neg).SL; % fz700
x_raw{3}{1}(:) = intersect_table_data(Data_900, SA_0).SL; % fz900
x_raw{3}{2}(:) = intersect_table_data(Data_900, SA_3neg).SL; % fz900
x_raw{3}{3}(:) = intersect_table_data(Data_900, SA_6neg).SL; % fz900
x_raw{4}{1}(:) = intersect_table_data(Data_1120, SA_0).SL; % fz1120
x_raw{4}{2}(:) = intersect_table_data(Data_1120, SA_3neg).SL; % fz1120
x_raw{4}{3}(:) = intersect_table_data(Data_1120, SA_6neg).SL; % fz1120

FX_raw = cell(4,1);
FX_raw{1}{1}(:) = intersect_table_data(Data_220, SA_0).FX; % fz220
FX_raw{1}{2}(:) = intersect_table_data(Data_220, SA_3neg).FX; % fz220
FX_raw{1}{3}(:) = intersect_table_data(Data_220, SA_6neg).FX; % fz220
FX_raw{2}{1}(:) = intersect_table_data(Data_700, SA_0).FX; % fz700
FX_raw{2}{2}(:) = intersect_table_data(Data_700, SA_3neg).FX; % fz700
FX_raw{2}{3}(:) = intersect_table_data(Data_700, SA_6neg).FX; % fz700
FX_raw{3}{1}(:) = intersect_table_data(Data_900, SA_0).FX; % fz900
FX_raw{3}{2}(:) = intersect_table_data(Data_900, SA_3neg).FX; % fz900
FX_raw{3}{3}(:) = intersect_table_data(Data_900, SA_6neg).FX; % fz900
FX_raw{4}{1}(:) = intersect_table_data(Data_1120, SA_0).FX; % fz1120
FX_raw{4}{2}(:) = intersect_table_data(Data_1120, SA_3neg).FX; % fz1120
FX_raw{4}{3}(:) = intersect_table_data(Data_1120, SA_6neg).FX; % fz1120

kappa_var = -0.25:0.001:0.25;
ones_vec = ones(length(kappa_var), 1);
FX_fit = cell(4,1);
[FX_fit{1}{1}(:)] = MF96_FXcomb(kappa_var, 0*pi/180*ones_vec, 0*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FX_fit{1}{2}(:)] = MF96_FXcomb(kappa_var, -3*pi/180*ones_vec, 0*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FX_fit{1}{3}(:)] = MF96_FXcomb(kappa_var, -6*pi/180*ones_vec, 0*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FX_fit{2}{1}(:)] = MF96_FXcomb(kappa_var, 0*pi/180*ones_vec, 0*ones_vec, mean(FZ_700.FZ)*ones_vec, tyre_coeffs);
[FX_fit{2}{2}(:)] = MF96_FXcomb(kappa_var, -3*pi/180*ones_vec, 0*ones_vec, mean(FZ_700.FZ)*ones_vec, tyre_coeffs);
[FX_fit{2}{3}(:)] = MF96_FXcomb(kappa_var, -6*pi/180*ones_vec, 0*ones_vec, mean(FZ_700.FZ)*ones_vec, tyre_coeffs);
[FX_fit{3}{1}(:)] = MF96_FXcomb(kappa_var, 0*pi/180*ones_vec, 0*ones_vec, mean(FZ_900.FZ)*ones_vec, tyre_coeffs);
[FX_fit{3}{2}(:)] = MF96_FXcomb(kappa_var, -3*pi/180*ones_vec, 0*ones_vec, mean(FZ_900.FZ)*ones_vec, tyre_coeffs);
[FX_fit{3}{3}(:)] = MF96_FXcomb(kappa_var, -6*pi/180*ones_vec, 0*ones_vec, mean(FZ_900.FZ)*ones_vec, tyre_coeffs);
[FX_fit{4}{1}(:)] = MF96_FXcomb(kappa_var, 0*pi/180*ones_vec, 0*ones_vec, mean(FZ_1120.FZ)*ones_vec, tyre_coeffs);
[FX_fit{4}{2}(:)] = MF96_FXcomb(kappa_var, -3*pi/180*ones_vec, 0*ones_vec, mean(FZ_1120.FZ)*ones_vec, tyre_coeffs);
[FX_fit{4}{3}(:)] = MF96_FXcomb(kappa_var, -6*pi/180*ones_vec, 0*ones_vec, mean(FZ_1120.FZ)*ones_vec, tyre_coeffs);

%%
legends_name = [];
plot_fitted_data_struct_combined(x_raw, FX_raw, kappa_var, FX_fit, '$\kappa [-]$', '$F_{X}$ [N]', ...
data_label, leg_angle, legends_name,'combined_longitudinal' , 'Combined longitudinal, $\gamma = 0 [deg]$', line_width, font_size_title, colors_vect)
%%

for i=1:length(kappa_var)
  [Gxa0(i), ~, ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), 0, 0, 0, tyre_coeffs);
  [Gxa3(i), ~, ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -3*pi/180, 0, 0, tyre_coeffs);
	[Gxa6(i), ~, ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -6*pi/180, 0, 0, tyre_coeffs);
  [Gxa8(i), ~, ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -8*pi/180, 0, 0, tyre_coeffs);
end

fig_Gxa = figure('Color', 'w');
plot(kappa_var, Gxa0, 'LineWidth', line_width);
hold on
plot(kappa_var, Gxa3, 'LineWidth', line_width);
plot(kappa_var, Gxa6, 'LineWidth', line_width);
plot(kappa_var, Gxa8, 'LineWidth', line_width);
hold off
grid on
xlabel('$\kappa [-]$', 'Interpreter', 'latex', 'FontSize', font_size)
ylabel('$G_{xa} [-]$', 'Interpreter', 'latex', 'FontSize', font_size)
title('Weighting function','Interpreter', 'latex', 'FontSize', font_size_title)
legend('$\alpha = 0 [deg]$', '$\alpha = -3 [deg]$', '$\alpha = -6 [deg]$', ...
'$\alpha = -8 [deg]$','Interpreter', 'latex', 'FontSize', font_size, ...
'location', 'southeast')
export_fig(fig_Gxa, 'images\fig_Gxa.png')

%% COMPARISON BETWEEEN THE TWO MEHTODS

ones_vec = ones(length(kappa_var),1);
sigma_x = kappa_var./(kappa_var + 1);
sigma_y0 = 0./(kappa_var + 1);
sigma_y3neg = -3*pi/180./(kappa_var + 1);
sigma_y6neg = -6*pi/180./(kappa_var + 1);
coeff0 = sigma_x./sqrt(sigma_x.^2 + sigma_y0.^2);
coeff3neg = sigma_x./sqrt(sigma_x.^2 + sigma_y3neg.^2);
coeff6neg = sigma_x./sqrt(sigma_x.^2 + sigma_y6neg.^2);

Fx = cell(4, 1);
Fx{1}{1}(:) = abs(coeff0).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 220*ones_vec,  tyre_coeffs);
Fx{1}{2}(:) = abs(coeff3neg).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 220*ones_vec,  tyre_coeffs);
Fx{1}{3}(:) = abs(coeff6neg).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 220*ones_vec,  tyre_coeffs);
Fx{2}{1}(:) = abs(coeff0).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 700*ones_vec,  tyre_coeffs);
Fx{2}{2}(:) = abs(coeff3neg).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 700*ones_vec,  tyre_coeffs);
Fx{2}{3}(:) = abs(coeff6neg).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 700*ones_vec,  tyre_coeffs);
Fx{3}{1}(:) = abs(coeff0).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 900*ones_vec,  tyre_coeffs);
Fx{3}{2}(:) = abs(coeff3neg).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 900*ones_vec,  tyre_coeffs);
Fx{3}{3}(:) = abs(coeff6neg).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 900*ones_vec,  tyre_coeffs);
Fx{4}{1}(:) = abs(coeff0).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 1120*ones_vec, tyre_coeffs);
Fx{4}{2}(:) = abs(coeff3neg).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 1120*ones_vec, tyre_coeffs);
Fx{4}{3}(:) = abs(coeff6neg).*MF96_FX0_vec(kappa_var, 0*ones_vec, 0*ones_vec, 1120*ones_vec, tyre_coeffs);
%%
% plot_fitted_data_struct_combined_sigma(x_raw, y_raw, x_fit, y_fit, x_fit_sig, y_fit_sig ...
%   label_x, label_y, data_label, leg_angle, name, plot_title, line_width, ...
%   font_size_title, colors_vector)
legend_names = [];
plot_fitted_data_struct_combined_sigma({}, {}, kappa_var, FX_fit, kappa_var, Fx, ...
'$\kappa [-]$', '$F_{x}$ [N]', data_label,leg_angle, legend_names,'combined_longitudinal_sigma' , ...
'Combined longitudinal, $\gamma =0$ [deg]', line_width, font_size_title, colors_vect)

PLON = {P1};
lbLON = {lb1};
ubLON = {ub1};
nameLON = {name1};
table_P0("combined_longitudinal", nameLON, PLON, lbLON, ubLON)

%% COMBINED LATERAL

%% Nominal conditions camber = 0, FZ = 220

TDataNominal = intersect_table_data(FZ_220, GAMMA_0);
FY_vec = TDataNominal.FY;
KAPPA_vec = TDataNominal.SL;
ALPHA_vec = TDataNominal.SA;
FZ_vec = TDataNominal.FZ; % 220 N
ones_vec = ones(length(FY_vec), 1);

clc
% rVy1, rVy4, rVy5, rVy6, rHy1, rBy1, rBy2, rBy3, rCy1
% [ -0.23 3.76 -0.09 28.37 0.02 14.16 13.29 -0.49 0.974];

P1 = 0.4*[1 1 1 1 1 1 1 1 1]; 
lb1 = [];
ub1 = [];
name1 = ["rVy1", "rVy4", "rVy5", "rVy6", "rHy1", "rBy1", "rBy2", "rBy3", "rCy1"];

[P_min,res_Fy,exitflag] = fmincon(@(P)resid_lateral_pure(P, FY_vec, KAPPA_vec, ALPHA_vec, 0 * ones_vec, FZ_vec, tyre_coeffs),...
P1,[],[],[],[],lb1,ub1);
P_min
res_Fy
RMS_Fy = sqrt(res_Fy*sum(FY_vec.^2)/length(FY_vec));
R2_Fy = 1 - res_Fy*sum(FY_vec.^2)/sum((FY_vec - mean(FY_vec)).^2);

tyre_coeffs.rVy1 = P_min(1);
tyre_coeffs.rVy4 = P_min(2);
tyre_coeffs.rVy5 = P_min(3);
tyre_coeffs.rVy6 = P_min(4);
tyre_coeffs.rHy1 = P_min(5);
tyre_coeffs.rBy1 = P_min(6);
tyre_coeffs.rBy2 = P_min(7);
tyre_coeffs.rBy3 = P_min(8);
tyre_coeffs.rCy1 = P_min(9);

TDataNominal_SA0 = intersect_table_data(SA_0, FZ_220, GAMMA_0);
TDataNominal_SA3neg = intersect_table_data(SA_3neg, FZ_220, GAMMA_0);
TDataNominal_SA6neg = intersect_table_data(SA_6neg, FZ_220, GAMMA_0);

KAPPA_vec = -0.25:0.001:0.25;
len = length(KAPPA_vec);
ones_vec = ones(len, 1);
ALPHA_vec0 = 0*ones_vec;
ALPHA_vec3neg = -3*pi/180*ones_vec;
ALPHA_vec6neg = -6*pi/180*ones_vec;
KAPPA_vec_struct = {KAPPA_vec, KAPPA_vec, KAPPA_vec};

FY_fit = cell(3, 1);
FY_fit{1} = MF96_FYcomb(KAPPA_vec, ALPHA_vec0, 0*ones_vec, mean(FZ_vec)*ones_vec, tyre_coeffs);
FY_fit{2} = MF96_FYcomb(KAPPA_vec, ALPHA_vec3neg, 0*ones_vec, mean(FZ_vec)*ones_vec, tyre_coeffs);
FY_fit{3} = MF96_FYcomb(KAPPA_vec, ALPHA_vec6neg, 0*ones_vec, mean(FZ_vec)*ones_vec, tyre_coeffs);

x_raw = {TDataNominal_SA0.SL, TDataNominal_SA3neg.SL, TDataNominal_SA6neg.SL};
y_raw = {TDataNominal_SA0.FY, TDataNominal_SA3neg.FY, TDataNominal_SA6neg.FY};

plot_label = {'$\alpha = 0 [deg]$', '$\alpha = -3 [deg]$', '$\alpha = -6 [deg]$'};
plot_fitted_data_struct(x_raw, y_raw, KAPPA_vec_struct, FY_fit, ...
'$\kappa$ [-]', '$F_{y}$ [N]', plot_label, 'fig_fit_nominal_conditions_FY_comb', ...
'Fitting in nominal conditions - Plots for different $\alpha$', line_width, ...
font_size_title, colors_vect);
  
  
%% Variable load
TDataFZ = GAMMA_0;
FY_vec = TDataFZ.FY;
KAPPA_vec = TDataFZ.SL;
ALPHA_vec = TDataFZ.SA;
FZ_vec = TDataFZ.FZ; % all FZ
ones_vec = ones(length(FY_vec), 1);

% rVy2 
P2 = [2]; 
lb2 = [];
ub2 = [];
name2 = ["rVy2"];

[P_min,res_Fy_varFz,exitflag] = fmincon(@(P)resid_lateral_varFz(P, ...
  FY_vec, KAPPA_vec, ALPHA_vec, 0 * ones_vec, FZ_vec, tyre_coeffs),...
P2,[],[],[],[],lb2,ub2);
P_min
res_Fy_varFz
RMS_Fy_varFz = sqrt(res_Fy_varFz*sum(FY_vec.^2)/length(FY_vec));
R2_Fy_varFz = 1 - res_Fy_varFz*sum(FY_vec.^2)/sum((FY_vec - mean(FY_vec)).^2);

tyre_coeffs.rVy2 = P_min(1);

TDataFZ_SA0 = intersect_table_data(SA_0, FZ_220, GAMMA_0);
TDataFZ_SA3neg = intersect_table_data(SA_3neg, FZ_220, GAMMA_0);
TDataFZ_SA6neg = intersect_table_data(SA_6neg, FZ_220, GAMMA_0);

KAPPA_vec = -0.25:0.001:0.25;
len = length(KAPPA_vec);
ones_vec = ones(len, 1);
ALPHA_vec0 = 0*ones_vec;
ALPHA_vec3neg = -3*pi/180*ones_vec;
ALPHA_vec6neg = -6*pi/180*ones_vec;
KAPPA_vec_struct = {KAPPA_vec, KAPPA_vec, KAPPA_vec};

FY_fit = cell(3, 1);
FY_fit{1} = MF96_FYcomb(KAPPA_vec, ALPHA_vec0, 0*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
FY_fit{2} = MF96_FYcomb(KAPPA_vec, ALPHA_vec3neg, 0*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
FY_fit{3} = MF96_FYcomb(KAPPA_vec, ALPHA_vec6neg, 0*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);

x_raw = {TDataFZ_SA0.SL, TDataFZ_SA3neg.SL, TDataFZ_SA6neg.SL};
y_raw = {TDataFZ_SA0.FY, TDataFZ_SA3neg.FY, TDataFZ_SA6neg.FY};
%%
plot_label = {'$\alpha = 0 [deg]$', '$\alpha = -3 [deg]$', '$\alpha = -6 [deg]$'};
plot_fitted_data_struct(x_raw, y_raw, KAPPA_vec_struct, FY_fit, ...
  '$\kappa$ [-]', '$F_{y}$ [N]', plot_label, 'fig_fit_variable_load_FY_comb', ...
  'Fitting with variable loads - Plots for different $\alpha$ and FZ=220[N]', line_width, ...
  font_size_title, colors_vect);
  
  %% Variable load and camber
TDataGamma = FZ_220;
FY_vec = TDataGamma.FY;
KAPPA_vec = TDataGamma.SL;
ALPHA_vec = TDataGamma.SA;
FZ_vec = TDataGamma.FZ; % Fz 220
GAMMA_vec = TDataGamma.IA; 
ones_vec = ones(length(FY_vec), 1);

% rVy3 
P3 = [2]; 
lb3 = [];
ub3 = [];
name3 = ["rVy3"];

[P_min,res_Fy_varGamma,exitflag] = fmincon(@(P)resid_lateral_varGamma(P, ...
  FY_vec, KAPPA_vec, ALPHA_vec, GAMMA_vec, FZ_vec, tyre_coeffs),...
P3,[],[],[],[],lb3,ub3);
P_min
res_Fy_varGamma
RMS_Fy_varGamma = sqrt(res_Fy_varGamma*sum(FY_vec.^2)/length(FY_vec));
R2_Fy_varGamma = 1 - res_Fy_varGamma*sum(FY_vec.^2)/sum((FY_vec - mean(FY_vec)).^2);

tyre_coeffs.rVy3 = P_min(1);

TDataFZ_SA0 = intersect_table_data(SA_0, FZ_220, GAMMA_0);
TDataFZ_SA3neg = intersect_table_data(SA_3neg, FZ_220, GAMMA_0);
TDataFZ_SA6neg = intersect_table_data(SA_6neg, FZ_220, GAMMA_0);

kappa_var = -0.25:0.01:0.25;

x_raw = cell(3, 1);
x_raw{1}{1}(:) = intersect_table_data(SA_0, FZ_220, GAMMA_0).SL;
x_raw{1}{2}(:) = intersect_table_data(SA_3neg, FZ_220, GAMMA_0).SL;
x_raw{1}{3}(:) = intersect_table_data(SA_6neg, FZ_220, GAMMA_0).SL;
x_raw{2}{1}(:) = intersect_table_data(SA_0, FZ_220, GAMMA_2).SL;
x_raw{2}{2}(:) = intersect_table_data(SA_3neg, FZ_220, GAMMA_2).SL;
x_raw{2}{3}(:) = intersect_table_data(SA_6neg, FZ_220, GAMMA_2).SL;
x_raw{3}{1}(:) = intersect_table_data(SA_0, FZ_220, GAMMA_4).SL;
x_raw{3}{2}(:) = intersect_table_data(SA_3neg, FZ_220, GAMMA_4).SL;
x_raw{3}{3}(:) = intersect_table_data(SA_6neg, FZ_220, GAMMA_4).SL;


y_raw = cell(3, 1);
y_raw{1}{1}(:) = intersect_table_data(SA_0, FZ_220, GAMMA_0).FY;
y_raw{1}{2}(:) = intersect_table_data(SA_3neg, FZ_220, GAMMA_0).FY;
y_raw{1}{3}(:) = intersect_table_data(SA_6neg, FZ_220, GAMMA_0).FY;
y_raw{2}{1}(:) = intersect_table_data(SA_0, FZ_220, GAMMA_2).FY;
y_raw{2}{2}(:) = intersect_table_data(SA_3neg, FZ_220, GAMMA_2).FY;
y_raw{2}{3}(:) = intersect_table_data(SA_6neg, FZ_220, GAMMA_2).FY;
y_raw{3}{1}(:) = intersect_table_data(SA_0, FZ_220, GAMMA_4).FY;
y_raw{3}{2}(:) = intersect_table_data(SA_3neg, FZ_220, GAMMA_4).FY;
y_raw{3}{3}(:) = intersect_table_data(SA_6neg, FZ_220, GAMMA_4).FY;

FY_fit = cell(4,1);
[FY_fit{1}{1}(:)] = MF96_FYcomb(kappa_var, 0*pi/180*ones_vec, mean(GAMMA_0.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FY_fit{1}{2}(:)] = MF96_FYcomb(kappa_var, -3*pi/180*ones_vec, mean(GAMMA_0.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FY_fit{1}{3}(:)] = MF96_FYcomb(kappa_var, -6*pi/180*ones_vec, mean(GAMMA_0.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FY_fit{2}{1}(:)] = MF96_FYcomb(kappa_var, 0*pi/180*ones_vec, mean(GAMMA_2.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FY_fit{2}{2}(:)] = MF96_FYcomb(kappa_var, -3*pi/180*ones_vec, mean(GAMMA_2.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FY_fit{2}{3}(:)] = MF96_FYcomb(kappa_var, -6*pi/180*ones_vec, mean(GAMMA_2.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FY_fit{3}{1}(:)] = MF96_FYcomb(kappa_var, 0*pi/180*ones_vec, mean(GAMMA_4.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FY_fit{3}{2}(:)] = MF96_FYcomb(kappa_var, -3*pi/180*ones_vec, mean(GAMMA_4.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);
[FY_fit{3}{3}(:)] = MF96_FYcomb(kappa_var, -6*pi/180*ones_vec, mean(GAMMA_4.IA)*ones_vec, mean(FZ_220.FZ)*ones_vec, tyre_coeffs);


data_label = ["$\gamma = 0 [deg]$", "$\gamma = 2 [deg]$", "$\gamma = 4 [deg]$"];
legends_name = ["Fit $\alpha = 0[deg]$", "Raw $\alpha$ = 0[deg]",...
"Fit $\alpha$ = -3[deg]", "Raw $\alpha$ = -3[deg]",...
"Fit $\alpha$ = -6[deg]", "Raw $\alpha$ = -6[deg]"];
plot_fitted_data_struct_combined(x_raw, y_raw, kappa_var, FY_fit, '$\kappa [-]$', '$F_{y}$ [N]', ...
data_label, leg_angle, legends_name,'combined_lateral' , 'Combined lateral, $FZ = 220 [N]$', line_width, font_size_title, colors_vect)


%% WEIGHTING FUNCTION
for i=1:length(kappa_var)
  [~, Gyk0(i), ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), 0, 0, 0, tyre_coeffs);
  [~, Gyk3(i), ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -3*pi/180, 0, 0, tyre_coeffs);
	[~, Gyk6(i), ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -6*pi/180, 0, 0, tyre_coeffs);
  [~, Gyk8(i), ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -8*pi/180, 0, 0, tyre_coeffs);
end

fig_Gyk = figure('Color', 'w');
plot(kappa_var, Gyk0, 'LineWidth', line_width);
hold on
plot(kappa_var, Gyk3, 'LineWidth', line_width);
plot(kappa_var, Gyk6, 'LineWidth', line_width);
plot(kappa_var, Gyk8, 'LineWidth', line_width);
hold off
grid on
xlabel('$\kappa [-]$', 'Interpreter', 'latex', 'FontSize', font_size)
ylabel('$G_{yk} [-]$', 'Interpreter', 'latex', 'FontSize', font_size)
title('Weighting function','Interpreter', 'latex', 'FontSize', font_size_title)
legend('$\alpha = 0 [deg]$', '$\alpha = -3 [deg]$', '$\alpha = -6 [deg]$', ...
'$\alpha = -8 [deg]$','Interpreter', 'latex', 'FontSize', font_size, ...
'location', 'northeast')
export_fig(fig_Gyk, 'images/fig_Gyk.png');
  
%% COMPARISON BETWEEEN THE TWO MEHTODS
ones_vec = ones(length(kappa_var),1);
sigma_x = kappa_var./(kappa_var + 1);
sigma_y0 = 0./(kappa_var + 1);
sigma_y3neg = -(3*pi/180)./(kappa_var + 1);
sigma_y6neg = -(6*pi/180)./(kappa_var + 1);
coeff0 = sigma_y0./sqrt(sigma_x.^2 + sigma_y0.^2);
coeff3neg = sigma_y3neg./sqrt(sigma_x.^2 + sigma_y3neg.^2);
coeff6neg = sigma_y6neg./sqrt(sigma_x.^2 + sigma_y6neg.^2);

Fy = cell(3, 1);
Fy{1}{1}(:) = abs(coeff0').*MF96_FY0_vec(0*ones_vec, 0*ones_vec, 0*ones_vec, 220*ones_vec, tyre_coeffs);
Fy{1}{2}(:) = abs(coeff3neg').*MF96_FY0_vec(0*ones_vec, -3*pi/180*ones_vec, 0*ones_vec, 220*ones_vec, tyre_coeffs);
Fy{1}{3}(:) = abs(coeff6neg').*MF96_FY0_vec(0*ones_vec, -6*pi/180*ones_vec, 0*ones_vec, 220*ones_vec, tyre_coeffs);
Fy{2}{1}(:) = abs(coeff0').*MF96_FY0_vec(0*ones_vec, 0*ones_vec, 2*pi/180*ones_vec, 220*ones_vec, tyre_coeffs);
Fy{2}{2}(:) = abs(coeff3neg').*MF96_FY0_vec(0*ones_vec, -3*pi/180*ones_vec, 2*pi/180*ones_vec, 220*ones_vec, tyre_coeffs);
Fy{2}{3}(:) = abs(coeff6neg').*MF96_FY0_vec(0*ones_vec, -6*pi/180*ones_vec, 2*pi/180*ones_vec, 220*ones_vec, tyre_coeffs);
Fy{3}{1}(:) = abs(coeff0').*MF96_FY0_vec(0*ones_vec, 0*ones_vec, 4*pi/180*ones_vec, 220*ones_vec, tyre_coeffs);
Fy{3}{2}(:) = abs(coeff3neg').*MF96_FY0_vec(0*ones_vec, -3*pi/180*ones_vec, 4*pi/180*ones_vec, 220*ones_vec, tyre_coeffs);
Fy{3}{3}(:) = abs(coeff6neg').*MF96_FY0_vec(0*ones_vec, -6*pi/180*ones_vec, 4*pi/180*ones_vec, 220*ones_vec, tyre_coeffs);
%%
% plot_fitted_data_struct_combined_sigma(x_raw, y_raw, x_fit, y_fit, x_fit_sig, y_fit_sig ...
%   label_x, label_y, data_label, leg_angle, name, plot_title, line_width, ...
%   font_size_title, colors_vector)
legend_names = ["Fit $\alpha = 0[deg]$","App. $\alpha = 0[deg]$", "Raw $\alpha$ = 0[deg]",...
"Fit $\alpha$ = -3[deg]","App. $\alpha$ = -3[deg]", "Raw $\alpha$ = -3[deg]",...
"Fit $\alpha$ = -6[deg]","App. $\alpha$ = -6[deg]", "Raw $\alpha$ = -6[deg]"];
plot_fitted_data_struct_combined_sigma(x_raw, y_raw, kappa_var, FY_fit, kappa_var, Fy, ...
'$\kappa [-]$', '$F_{y}$ [N]', data_label,leg_angle, legend_names,'combined_lateral_sigma' , ...
'Combined lateral, $FZ = 220 [N]$', line_width, font_size_title, colors_vect)


fprintf("--- Combined Fx ---\n")
fprintf("Nominal conditions: %6.3f\n", res_Fx);
fprintf("Nominal conditions: %6.3f\n", R2_Fx);
fprintf("Nominal conditions: %6.3f\n", RMS_Fx);

fprintf("--- Combined Fy ---\n")
fprintf("Nominal conditions: %6.3f\n", res_Fy);
fprintf("Variable load: %6.3f\n", res_Fy_varFz);
fprintf("Variable camber: %6.3f\n", res_Fy_varGamma);
fprintf('\n')
fprintf("Nominal conditions: %6.3f\n", R2_Fy);
fprintf("Variable load: %6.3f\n", R2_Fy_varFz);
fprintf("Variable camber: %6.3f\n", R2_Fy_varGamma);
fprintf('\n')
fprintf("Nominal conditions: %6.3f\n", RMS_Fy);
fprintf("Variable load: %6.3f\n", RMS_Fy_varFz);
fprintf("Variable camber: %6.3f\n", RMS_Fy_varGamma);


%%
PLAT = {P1 P2 P3};
lbLAT = {lb1 lb2 lb3};
ubLAT = {ub1 ub2 ub3}; 
nameLAT = {name1 name2 name3};
table_P0("combined_lateral", nameLAT, PLAT, lbLAT, ubLAT)
%% Save tyre data structure to mat file
name = ["resFx", "RtwoFx", "RMSFx", "resFy", "resFyvarFz", "resFyvarGamma", "RtwoFy", "RtwoFyvarFz", "RtwoFyvarGamma", "RMSFy", "RMSFyvarFz", "RMSFyvarGamma"];

data = [res_Fx, R2_Fx, RMS_Fx, res_Fy, res_Fy_varFz, res_Fy_varGamma, R2_Fy, R2_Fy_varFz, R2_Fy_varGamma, RMS_Fy, RMS_Fy_varFz, RMS_Fy_varGamma];


write_latex_macro('results_to_lates_Fcomb.tex', name, data, 'w');
save(['tyre_' struct_name,'.mat'],'tyre_coeffs');