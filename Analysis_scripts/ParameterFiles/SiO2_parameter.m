% Parameter file for Al/Interface/SiO2/Interface/Si

function [SysParam] = SiO2_parameter() %Apr19Fits_Si_025nmSiO2_2MHz() %This is a test file

%% DATA FILE NAME

SysParam.filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a11_Si_SiO2_50nm_NDto800nmAdded_corrected_150623_200157_FIN.mat';

%% SAMPLE PROPERTIES (starting with top layer )
% Al / Interface / SiO2 / Interface / Si:
% Al and Si properties are the same as for "Apr19Fits_Si_PECVD_noSiO2_2MHz" file.

SysParam.Lambda  = [237 0.126 0.781 999 124];  % Thermal conductivities (W m^-1 K^-1)
SysParam.C       = [2.49 2.42 0.949 2.42 1.63]*1e6;  % Volumetric heat capacities (J m^-3 K^-1)
SysParam.h       = [80 1 56.2 1 1e6]*1e-9;  % Thicknesses (m)  

SysParam.eta     = ones(1,numel(SysParam.Lambda)); % Anisotropy parameter eta=kx/ky;
SysParam.X_heat  = ((1:5:36)*1e-9)';  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
SysParam.X_temp  = 1e-9;   % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
SysParam.AbsProf = exp(-4*pi*4.9*SysParam.X_heat/400e-9);  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized

%% SETUP PROPERTIES
SysParam.f        = 2.0e6;  % Laser modulation frequency (Hz)
SysParam.r_pump   = 0.5*72.0e-6;  % Pump 1/e^2 radius (m) % previously 0.5*72.0e-6;
SysParam.r_probe  = 0.5*23.5e-6;  % Probe 1/e^2 radius, (m) % previously 0.5*29.4e-6;
SysParam.tau_rep  = 1/(75.8e6);   % Laser repetition period (s) 
SysParam.P_pump   = 0.8*21.0e-3; % 0.99*0.3*37e-3;   % absorbed pump power (transmission of objective * absorbance * pump power) 
SysParam.P_probe  = 0.4*50.0e-3; % 0.8*0.3*17.5e-3;   % absorbed probe power (transmission of objective * absorbance * probe power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
SysParam.tdelay_model = logspace(log10(10e-12),log10(4e-9),60)'; % time delays for model curve (s)
% VKORN :: changed to 60 from 30 in logspace().

%% FIT SPECIFICATIONS
% Choose layer indices; respective parameters are then adjusted in the fit;
SysParam.FITNLambda = [2 3]; % [2 3 4]
SysParam.FITNC = [3]; 
SysParam.FITNh = [];
% Choose range of time delays to fit (s)
SysParam.tdelay_min = 120e-12; % before Apr19: 30e-12; 
SysParam.tdelay_max = 5000e-12; %
