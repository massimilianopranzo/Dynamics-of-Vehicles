
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
font_size = 18;       % font size for the labels
font_size = 15; % font size for the titles
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

colors_vect = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; ...
               [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; ...
               [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; ...
               [0.6350 0.0780 0.1840]];

