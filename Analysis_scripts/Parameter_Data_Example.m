function [SysParam] = Parameter_Data_Example()

SysParam.filename = 'Data_Example.mat';

%% SAMPLE PROPERTIES (starting with top layer )
SysParam.Lambda  = [45  0.3 40];  % Thermal conductivities (W m^-1 K^-1)
SysParam.C       = [2.9 0.1 3.36]*1e6;  % Volumetric heat capacities (J m^-3 K^-1)
SysParam.h       = [60 1 1e6]*1e-9;  % Thicknesses (m)   
SysParam.eta     = ones(1,numel(SysParam.Lambda)); % Anisotropy parameter eta=kx/ky;
SysParam.X_heat  = ((1:5:36)*1e-9)';  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
SysParam.X_temp  = 1e-9;   % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
SysParam.AbsProf = exp(-4*pi*4.9*SysParam.X_heat/785e-9);  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized

%% SETUP PROPERTIES
SysParam.f        = 10.1e6;  % Laser modulation frequency (Hz)
SysParam.r_pump   = 5.3e-6;  % Pump 1/e^2 radius (m)
SysParam.r_probe  = 5.3e-6;  % Probe 1/e^2 radius, (m)
SysParam.tau_rep  = 1/80e6;   % Laser repetition period (s) 
SysParam.P_pump   = 0.8*0.3*7e-3;   % absorbed pump power (transmission of objective * absorbance * pump power) 
SysParam.P_probe  = 0.8*0.3*6e-3;   % absorbed pump power (transmission of objective * absorbance * pump power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
SysParam.tdelay_model = logspace(log10(10e-12),log10(4e-9),30)'; % time delays for model curve (s)

%% FIT SPECIFICATIONS
% Choose layer indices; respective parameters are then adjusted in the fit;
SysParam.FITNLambda = [2,3];
SysParam.FITNC = []; 
SysParam.FITNh = [];
% Choose range of time delays to fit (s)
SysParam.tdelay_min = 100e-12; 
SysParam.tdelay_max = 4000e-12;
