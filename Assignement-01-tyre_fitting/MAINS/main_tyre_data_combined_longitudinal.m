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
plots_raw_data;

%% Sort data
sorting_data;
% figure 2
plot_sorted_data;

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
P0 = [10 10 2 9];
lb = [];
ub = [];

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
