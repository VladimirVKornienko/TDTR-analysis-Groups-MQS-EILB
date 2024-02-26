function [SysParam] = dummy_Au_SiO2_Si()

% This file is created in order to be used with "_SIM" routine, to
% find the best possible structure for recovering the Au/Ti transducer
% properties.

%% DATA FILE NAME

SysParam.filename = 'Z:\ForMATLAB\data\RefAuGlass1MHz@1MHz_011123_210559_FIN.mat';
% <<< does not match any 'real' data, just something; it will not be used...%

%% SAMPLE PROPERTIES (starting with top layer )
% Au / SiO2 / Si: no need to model the thermal interfaces, if there is SiO2
% which is already a poor heat conductor.

% E.g.: Au (50 nm) / SiO2 (0 / 50 / 100 nm) / Si (substr.)

SysParam.Labels = ["Au", "SiO2","Si (substr.)"]; % Labels for each element in the lists below. Only used in SIM plots

%% Initial guess for the parameters (table values):
% From the 2 MHz datas -- initial guess:
SysParam.Lambda  = [315 1.3 124];  % Thermal conductivities (W m^-1 K^-1)
SysParam.C       = [2.49 1.86 1.63]*1e6;  % Volumetric heat capacities (J m^-3 K^-1)
SysParam.h       = [50.0 5 1e6]*1e-9;  % Thicknesses (m)


%% other parameters
SysParam.eta     = ones(1,numel(SysParam.Lambda)); % Anisotropy parameter eta=kx/ky;
SysParam.X_heat  = ((1:5:36)*1e-9)';  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
SysParam.X_temp  = 1e-9;   % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
SysParam.AbsProf = exp(-4*pi*4.9*SysParam.X_heat/400e-9);  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized

%% SETUP PROPERTIES
SysParam.f        = 1.0e6;  % Laser modulation frequency (Hz)
SysParam.r_pump   = 0.5*56.1e-6;  % Pump 1/e^2 radius (m)


% temp Nov. 17 >>
% SysParam.r_pump   = 0.5*35.0e-6;  % Pump 1/e^2 radius (m)
% <<<


SysParam.r_probe  = 0.5*20.3e-6;  % Probe 1/e^2 radius, (m)
SysParam.tau_rep  = 1/(75.8e6);   % Laser repetition period (s) 
SysParam.P_pump   = 0.8*21.0e-3; % 0.99*0.3*37e-3;   % absorbed pump power (transmission of objective * absorbance * pump power) 
SysParam.P_probe  = 0.4*50.0e-3; % 0.8*0.3*17.5e-3;   % absorbed probe power (transmission of objective * absorbance * probe power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
%% 23.02.2024: >>>
%SysParam.tdelay_model = logspace(log10(10e-12),log10(4e-9),60)'; % time delays for model curve (s)
% VKORN :: changed to 60 from 30 in logspace().
SysParam.tdelay_model = logspace(log10(50e-12),log10(4e-9),30)'; % time delays for model curve (s)
%% <<<<<

%% FIT SPECIFICATIONS
% Choose layer indices; respective parameters are then adjusted in the fit;
SysParam.FITNLambda = [1 2 3]; % [2 3 4]
SysParam.FITNC = [1 2 3]; 
SysParam.FITNh = [1 2];

SysParam.FITNLambda = [1 2];
SysParam.FITNC = []; 
SysParam.FITNh = [];

% Choose range of time delays to fit (s)
% SysParam.tdelay_min = 75e-12; % before Apr19: 30e-12; 
SysParam.tdelay_min = 200e-12; % <<< upd. 24.11.2023 <<< ; 
SysParam.tdelay_max = 5000e-12;