% Parameter file for Al/Interface/CuI/Interface/[Glass or Si]

function [SysParam] = CuI_parameter() %Apr19Fits_Si_025nmSiO2_2MHz() %This is a test file

%% DATA FILE NAME

% CuI thicknesses: 0deg: 1010, 70deg: 290, 75deg: 262, 80deg: 137
projectPath = 'C:\Users\immonel1\Downloads\Analysis_scripts';
dataFolder = '\2023_07_12 CuI new\f03_preprocessed_data';
%SysParam.filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a04_CuI_point2_redone_150623_165645_FIN.mat';
%SysParam.filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a06_CuI_70degr_GLAD_point2_corrected_150623_175319_FIN.mat';
%SysParam.filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a08_CuI_75degr_glad_point2_corrected_150623_184625_FIN.mat';
%SysParam.filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a09_CuI_80degr_glad_point1_150623_191225_FIN.mat'; 
%filename = '\f03_CuI_0deg_point1_120723_163922_FIN.mat';
%filename = '\f04_CuI_0deg_point2_120723_172648_FIN.mat';
%filename = '\f05_CuI_80deg_1min_point1_120723_175747_FIN.mat';
filename = '\f06_CuI_80deg_1min_point2_120723_182727_FIN.mat';
%filename = '\f07_CuI_80deg_2min_point1_120723_185415_FIN.mat';
%filename = '\f08_CuI_80deg_2min_point2_120723_192033_FIN.mat';
SysParam.filename = [projectPath dataFolder filename];

transducer = "Au";

%% SAMPLE PROPERTIES (starting with top layer )
% Al / Interface / CuI / Interface / Glass:
% Al and Si properties are the same as for "Apr19Fits_Si_PECVD_noSiO2_2MHz" file.

if transducer == "Al"
    SysParam.Labels = ["Al", "Al/CuI", "CuI", "CuI/Si", "Si"]; % Labels for each element in the lists below. Only used in SIM plots
    SysParam.Lambda  = [237 .8 0.229 999 124];  % Thermal conductivities (W m^-1 K^-1) [237 0.129 1.68 999 1.14];
    SysParam.C       = [1.3 2.42 13.5 2.42 1.63]*1e6;  % Volumetric heat capacities (J m^-3 K^-1) [2.42 2.42 1.61 2.42 1.85]*1e6;
    SysParam.h       = [80 1 66 1 1e6]*1e-9;  % Thicknesses (m)  [78.25 1 137 1 1e6]*1e-9;
elseif transducer == "Au"
    SysParam.Labels = ["Au", "Ti", "CuI", "CuI/Si", "Si"]; % Labels for each element in the lists below. Only used in SIM plots
    SysParam.Lambda  = [320 17 1.68 999 124];  % Thermal conductivities (W m^-1 K^-1) [237 0.129 1.68 999 1.14];
    SysParam.C       = [2.49 2.64 1.61 2.42 1.63]*1e6;  % Volumetric heat capacities (J m^-3 K^-1) [2.42 2.42 1.61 2.42 1.85]*1e6;
    SysParam.h       = [150 5 66 1 1e6]*1e-9;  % Thicknesses (m)  [78.25 1 137 1 1e6]*1e-9;
else
    disp("ERROR! CHECK VALUE OF VARIABLE 'transducer'")
end

SysParam.eta     = ones(1,numel(SysParam.Lambda)); % Anisotropy parameter eta=kx/ky;
SysParam.X_heat  = ((1:5:36)*1e-9)';  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
SysParam.X_temp  = 1e-9;   % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
SysParam.AbsProf = exp(-4*pi*4.9*SysParam.X_heat/400e-9);  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized

% testing what happens if AbsProf = [1,0,0,0, ...]
%SysParam.AbsProf = zeros(size(SysParam.X_heat));
%5ysParam.AbsProf(1) = 1;

%% SETUP PROPERTIES
SysParam.f        = 2.0e6;  % Laser modulation frequency (Hz)
SysParam.r_pump   = 0.5*95.0e-6;  % Pump 1/e^2 radius (m) % previously 0.5*72.0e-6;
SysParam.r_probe  = 0.5*15.0e-6;  % Probe 1/e^2 radius, (m) % previously 0.5*29.4e-6;
SysParam.tau_rep  = 1/(75.8e6);   % Laser repetition period (s) 
SysParam.P_pump   = 0.8*21.0e-3; % 0.99*0.3*37e-3;   % absorbed pump power (transmission of objective * absorbance * pump power) 
SysParam.P_probe  = 0.4*50.0e-3; % 0.8*0.3*17.5e-3;   % absorbed probe power (transmission of objective * absorbance * probe power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
SysParam.tdelay_model = logspace(log10(10e-12),log10(4e-9),60)'; % time delays for model curve (s)
% VKORN :: changed to 60 from 30 in logspace().

%% FIT SPECIFICATIONS
% Choose layer indices; respective parameters are then adjusted in the fit;
SysParam.FITNLambda = [1 2 3]; % [2 3 4]
SysParam.FITNC = [1 2 3]; 
SysParam.FITNh = [1 2 3];
% Choose range of time delays to fit (s)
SysParam.tdelay_min = 60e-12; % before Apr19: 30e-12; 120e-12;
SysParam.tdelay_max = 2000e-12;%2e-9;%
