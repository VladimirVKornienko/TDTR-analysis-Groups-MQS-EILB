function [SysParam] = curr_Nov07_5MHz_CuI_70deg()

%% DATA FILE NAME

SysParam.filename = 'Z:\ForMATLAB\data\Cu_70glad_30sec_5MHz_NoOsc@5MHz_011123_225301_FIN.mat';

%% SAMPLE PROPERTIES (starting with top layer )
% CuI 80 deg. GLAD (116) / interface (1 nm) / Au (50 nm) / Ti (5 nm) / SiO2 (substr.)
% SysParam.Labels = ["Au", "Ti", "CuI", "CuI/SiO2", "SiO2 (substr.)"]; % Labels for each element in the lists below. Only used in SIM plots

% update 24.02.2024:
% maybe, we can drop the thermal interfaces for now, let's try just the
% CuI:
% SysParam.Labels = ["Au", "CuI", "SiO2 (substr.)"]; % Labels for each element in the lists below. Only used in SIM plots
%
% Neah. Seems like we still need the Au/CuI interface... 
% Still no. Let's try the interface between the CuI and the SiO2:
% SysParam.Labels = ["Au", "CuI", "CuI/SiO2", "SiO2 (substr.)"]; % Labels for each element in the lists below. Only used in SIM plots
% ok, this gives better results. For 1 nm, we get lambda = 0.002 for the
% interface.
%
% ok. seems, we need the 2 interfaces still...
% SysParam.Labels = ["Au", "Au/CuI", "CuI", "CuI/SiO2", "SiO2 (substr.)"]; % Labels for each element in the lists below. Only used in SIM plots

% after fitting: CuI / SiO2 interface is not needed %
SysParam.Labels = ["Au", "Au/CuI", "CuI", "SiO2 (substr.)"]; % Labels for each element in the lists below. Only used in SIM plots

%% Initial guess for the parameters:
% SysParam.Lambda  = [315 0.08 1.68 0.02 1.24];  % Thermal conductivities (W m^-1 K^-1)
% SysParam.C       = [2.49 2.49 1.55 2.49 1.86]*1e6;  % Volumetric heat capacities (J m^-3 K^-1)
% SysParam.h       = [54 1 196 1 1e6]*1e-9;  % Thicknesses (m)

SysParam.Lambda  = [315 0.08 3.0 1.24];  % Thermal conductivities (W m^-1 K^-1)
SysParam.C       = [2.49 2.49 1.55 1.86]*1e6;  % Volumetric heat capacities (J m^-3 K^-1)
SysParam.h       = [56 1 196 1e6]*1e-9;  % Thicknesses (m)



%% Temp for fitting:
SysParam.Lambda  = [310.3 0.367 1.68 1.152];  % Thermal conductivities (W m^-1 K^-1)
SysParam.C       = [2.363 2.64 1.55 1.871]*1e6;  % Volumetric heat capacities (J m^-3 K^-1)
SysParam.h       = [50.1 5 196 1e6]*1e-9;  % Thicknesses (m)
% <<< %


%% upd. 26.02.2024: %%
% TOLERANCES / SPREADs: >>> %
SysParam.rangeLambda  = 1000*[30 100 100 2];  % Thermal conductivities (W m^-1 K^-1)
SysParam.rangeC       = 1000*[0.4 0.0001 100 0.4]*1e6;  % Volumetric heat capacities (J m^-3 K^-1)
SysParam.rangeh       = 100*[5.0 0.0001 10.0 0.0001]*1e-9;  % Thicknesses (m)
% <<< %


%% other parameters
SysParam.eta     = ones(1,numel(SysParam.Lambda)); % Anisotropy parameter eta=kx/ky;
SysParam.X_heat  = ((1:5:36)*1e-9)';  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
SysParam.X_temp  = 1e-9;   % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
SysParam.AbsProf = exp(-4*pi*4.9*SysParam.X_heat/400e-9);  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized

%% SETUP PROPERTIES
SysParam.f        = 5.0e6;  % Laser modulation frequency (Hz)
SysParam.r_pump   = 0.5*56.1e-6;  % Pump 1/e^2 radius (m)
SysParam.r_probe  = 0.5*20.3e-6;  % Probe 1/e^2 radius, (m)
SysParam.tau_rep  = 1/(75.8e6);   % Laser repetition period (s) 
SysParam.P_pump   = 0.8*21.0e-3; % 0.99*0.3*37e-3;   % absorbed pump power (transmission of objective * absorbance * pump power) 
SysParam.P_probe  = 0.4*50.0e-3; % 0.8*0.3*17.5e-3;   % absorbed probe power (transmission of objective * absorbance * probe power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
SysParam.tdelay_model = logspace(log10(50e-12),log10(4e-9),30)'; % time delays for model curve (s)
% VKORN :: changed to 60 from 30 in logspace().

%% FIT SPECIFICATIONS
% Choose layer indices; respective parameters are then adjusted in the fit;
SysParam.FITNLambda = [2 3 4];
SysParam.FITNC = [1 3 4];
SysParam.FITNh = [];

SysParam.FITNLambda = [2 3];
SysParam.FITNC = [3]; 
SysParam.FITNh = [];

% Choose range of time delays to fit (s)
SysParam.tdelay_min = 200e-12; % before Apr19: 30e-12; 
SysParam.tdelay_max = 5000e-12;
