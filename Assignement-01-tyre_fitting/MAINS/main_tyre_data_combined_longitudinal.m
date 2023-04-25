%% Initialisation
clc
clearvars 
close all 

addpath('utilities\')
addpath('MAINS\..')
addpath('tyre_lib\')

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
  tyre_coeffs.rHx1 = 0;
  tyre_coeffs.rBx1 = 0;
  tyre_coeffs.rCx1 = 0;
  tyre_coeffs.rBx2 = 0;
   
else
  error('Run the genertion of FX, FY and MZ parameters first!');
end



% rHx1, rBx1, rCx1 rBx2
% -27.1725    5.5548   11.0463   -1.4499
P0 = [-27 5.5 11 -1.5];
lb = [-30 4 10 -2];
ub = [-20 7 15 0];

FX_vec = tyre_data.FX;
KAPPA_vec = tyre_data.SL;
ALPHA_vec = tyre_data.SA;
GAMMA_vec = tyre_data.IA;
FZ_vec = tyre_data.FZ;

[P_min,fval,exitflag] = fmincon(@(P)resid_long(P, FX_vec, KAPPA_vec, ALPHA_vec, GAMMA_vec, FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

P_min
fval

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

plot_fitted_data_struct(x_raw, FX_raw, x_raw, FX_fit, '$\kappa [-]$', 'FX [N]', ...
  data_label, 'combined_longitudinal' , 'Combined longitudinal, $\gamma = 0 [deg]$', line_width, font_size_title, colors_vect)


kappa_var = -0.5:0.001:0.5;
for i=1:length(kappa_var)
  [Gxa0(i), ~, ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), 0, 0, 0, tyre_coeffs);
  [Gxa3(i), ~, ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -3*pi/180, 0, 0, tyre_coeffs);
	[Gxa6(i), ~, ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -6*pi/180, 0, 0, tyre_coeffs);
  [Gxa8(i), ~, ~] = MF96_FXFYCOMB_coeffs_eqns(kappa_var(i), -8*pi/180, 0, 0, tyre_coeffs);
end

figure('Color', 'w')
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