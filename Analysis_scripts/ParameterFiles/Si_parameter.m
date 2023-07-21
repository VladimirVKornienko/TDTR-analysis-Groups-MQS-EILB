% Parameter file for Al/Interface/Si

function [SysParam] = Si_parameter() %Apr19Fits_Si_025nmSiO2_2MHz() %This is a test file

%% DATA FILE NAME
projectPath = 'C:\Users\immonel1\Downloads\Analysis_scripts';
dataFolder = '\2023_07_12 CuI new\f03_preprocessed_data';
%SysParam.filename = 'dataForFits\SiO2_2MHz_180423_121206_FIN_NoNaNs.mat';
%SysParam.filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a01_Si_Al80nm_point1_sequence1_150623_131327_FIN.mat';
%SysParam.filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a02_Si_second_run_different_sequence_corrected_150623_141649_FIN.mat';
%SysParam.filename = '\Tarmo_data\2023_07_07\Si_2MHz_070723_162924.mat';
%filename = '\f01_Si_point1_120723_153555_FIN.mat';
filename = '\f02_Si_point2_120723_155842_FIN.mat';
SysParam.filename = [projectPath dataFolder filename];

%% SAMPLE PROPERTIES (starting with top layer )
% Al / Interface / Si:
% Al and Si properties are the same as for "Apr19Fits_Si_PECVD_noSiO2_2MHz" file.

SysParam.Labels = ["Al", "Al/Si", "Si"]; % Labels for each element in the lists below. Only used in SIM plots
SysParam.Lambda  = [237 0.129 124];  % Thermal conductivities (W m^-1 K^-1) for Al / Interface / Si [237 0.15 124];
SysParam.C       = [2.42 2.42 1.63]*1e6;  % Volumetric heat capacities (J m^-3 K^-1) %2.42  for Al / Interface / Si  [2.42 2.42 1.63]*1e6;
SysParam.h       = [80 1 1e6]*1e-9;  % Thicknesses (m)   for Al / Interface / Si [80 1 1e6]*1e-9;

SysParam.eta     = ones(1,numel(SysParam.Lambda)); % Anisotropy parameter eta=kx/ky;
SysParam.X_heat  = ((1:5:36)*1e-9)';  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
SysParam.X_temp  = 1e-9;   % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
SysParam.AbsProf = exp(-4*pi*4.9*SysParam.X_heat/400e-9);  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized

% testing what happens if AbsProf = [1,0,0,0, ...]
%SysParam.AbsProf = zeros(size(SysParam.X_heat));
%SysParam.AbsProf(1) = 1;

%% SETUP PROPERTIES
SysParam.f        = 2.0e6;  % Laser modulation frequency (Hz)
SysParam.r_pump   = 0.5*95.0e-6;  % Pump 1/e^2 radius (m) % previously 0.5*72.0e-6;
SysParam.r_probe  = 0.5*15e-6;  % Probe 1/e^2 radius, (m) % previously 0.5*29.4e-6;
SysParam.tau_rep  = 1/(75.8e6);   % Laser repetition period (s) 
SysParam.P_pump   = 0.8*21.0e-3; % 0.99*0.3*37e-3;   % absorbed pump power (transmission of objective * absorbance * pump power) 
SysParam.P_probe  = 0.4*50.0e-3; % 0.8*0.3*17.5e-3;   % absorbed probe power (transmission of objective * absorbance * probe power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
SysParam.tdelay_model = logspace(log10(10e-12),log10(4e-9),60)'; % time delays for model curve (s)
% VKORN :: changed to 60 from 30 in logspace().

%% FIT SPECIFICATIONS
% Choose layer indices; respective parameters are then adjusted in the fit;
% For auto fit: [2], [1] and []. For manual [1 2], [1], [1]
SysParam.FITNLambda = [2];%[1 2]; % [2 3 4]
SysParam.FITNC = [1];%1]; 
SysParam.FITNh = [];%1];
% Choose range of time delays to fit (s)
SysParam.tdelay_min = 120e-12; % 75e-12 for a1. 120e-12 for a2 % before Apr19: 30e-12; 
SysParam.tdelay_max = 5000e-12;
