% Quick upgrade to fit to the ratio of ratios between the two frequencies.
% 23.02.2024

% It calls the "regular" subroutine SUB_C, gets the vector of
% <ratio_model>, computes the R.M.S.E. of the between (R1/R2)_data and
% (R1/R2)_model. It then returns R.M.S.E. ("ZratioERR") to
% provide a value for fminsearch() to work on. 

% <<<<< 23.02.2024 <<<<< %


% Input parameters:

% In comparison to "SUB_C": >>>
% three additional arguments, that go at the
% start of the new function:
%
% f2_Ratio_data, f2_tdelay, f2_f.
%
% <<<< %


% figNum - for "figure(figNum)".

%   X           - New values of the fit parameters
%   ratio_data  - TDTR ratio data measured
%   tdelay      - time delay (ps)
%   rau_rep     - repetition rate of laser (s)
%   f           - modulation frequency (Hz)
%   Lambda      - vector of thermal conductivities, Lambda(1) = top surface,(W/m/K)
%   C           - vector of volumetric heat capacity (J/m^3/K)
%   h           - thicknesses of each layer (layer N will NOT be used, semiinfinite)
%   r_pump      -  pump 1/e^2 radius (m)
%   r_probe     - probe 1/e^2 radius (m)
%   P_pump      - absorbed pump laser power (W) ...doesn't effect ratio
%   nnodes      - number of nodes used for numerical integration
%   FITNLambda  - row vector of indices of layers of which the thermal conductivity is adjusted in the fit
%   FITNC       - row vector of indices of layers of which the heat capacity is adjusted in the fitting
%   FITNh       - row vector of indices of layers of which the layer thickness is adjusted in the fitting
%   X_heat      - temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
%   X_temp      - depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)

%---------------------------- BEGIN CODE ----------------------------------

function [ZratioERR] = TDTR_Bidirectional_SUB_VKorn_D(f2_Ratio_data,f2_tdelay,f2_f, figNum, X,Ratio_data,tdelay,tau_rep,f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf)

% "SUB_C" returns: % function [Z,Ratio_model] % 
[subD_Z1,file1_Ratio_model] = TDTR_Bidirectional_SUB_C(figNum, X,Ratio_data,tdelay,tau_rep,f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf);
%fprintf("%5.3f\t",file1_Ratio_model);
[subD_Z2,file2_Ratio_model] = TDTR_Bidirectional_SUB_C(figNum+1, X,f2_Ratio_data,f2_tdelay,tau_rep,f2_f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf);
%fprintf("%5.3f\t",file2_Ratio_model);

% in "SUB_C": Z = ((Ratio_model-Ratio_data)./Ratio_model).^2;
% here: (R1/R2)_model - (R1/R2)_data / (R1/R2)_model .^2;

subD_model = (file1_Ratio_model ./ file2_Ratio_model);
subD_data = (Ratio_data ./ f2_Ratio_data);

ZratioERR = (  (subD_model - subD_data) ./ subD_model  ).^2;
ZratioERR = sqrt(sum(ZratioERR))/length(ZratioERR);

% IMP! Add the original files' fit errors to the ratio error, %
% tune the coefficients, if needed... %
subD_c1empyric = 0.2;
subD_c2empyric = 0.2;
ZratioERR = sqrt( subD_c1empyric*(subD_Z1)^2 + subD_c2empyric*(subD_Z2)^2 + (ZratioERR)^2);
fprintf("Total error: %f\n",ZratioERR);
% <<< %

% fprintf("%f",ZratioERR);

figure(figNum+2)  % +0 is file1, +1 is file2, +2 is for this, ratio, graph
clf
%semilogx(tdelay, -real(Ts)./imag(Ts),'g')
hold on
% semilogx(tdelay,Ratio_data,'ob',tdelay,Ratio_model,'k');
semilogx(tdelay,subD_data,'ob',tdelay,subD_model,'k');
pause(0.01)

