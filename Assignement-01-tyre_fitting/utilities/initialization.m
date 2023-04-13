%% Initialisation
clc
clearvars 
close all 

% Choice of the dataset
data_set_path = 'dataset/';
data_set = 'Goodyear_B1464run58';  
% tyre geometric data:
% Goodyear D2704 20.0x7.0-13
% 20 diameter in inches
% 7.0 section width in inches
% tread width in inches
diameter = 20*2.56; % [cm]
Fz0 = 220;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***

% Set valeus for the plots
font_size_title = 20; % font size for the titles
line_width = 2;       % line width for the plots

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  16)
set(0,'DefaultLegendFontSize',16)


addpath('tyre_lib\')

to_rad = pi/180;
to_deg = 180/pi;

