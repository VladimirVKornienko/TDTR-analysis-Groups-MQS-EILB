% This function is part of the "Bidirectional" TDTR analysis.
% This function is used for fitting the TDTR_Surface model to the TDTR ratio measured: The main program tries to minimize "Z" by optimizing the variable(s) X.
% This function is based on TDTR_FIT_V4.m published by Joseph P. Feser under http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html (12-September-2012).
% This function can handle varying pump size over delay time thanks to Greg Hohensee.

% Input paameters:

% ADDED: figNum - for "figure(figNum)".

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

function [Z,Ratio_model] = TDTR_Bidirectional_SUB_C(figNum, X,Ratio_data,tdelay,tau_rep,f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf)

% assign fit variable X to material parameters 
    for i = 1:length(FITNLambda)
        Lambda(FITNLambda(i)) = X(i);   
    end
    for i = 1:length(FITNC)
        C(FITNC(i)) = X(length(FITNLambda)+i);   
    end
    for i = 1:length(FITNh)
        h(FITNh(i)) = X(length(FITNLambda)+length(FITNC)+i);
    end

[Ts,~] = TDTR_Bidirectional_SUB_B(tdelay,tau_rep,f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,X_heat,X_temp);

Tin_model = real(Ts)*AbsProf./(ones(size(AbsProf))'*AbsProf);
Tout_model = imag(Ts)*AbsProf./(ones(size(AbsProf))'*AbsProf);
Ratio_model = -Tin_model./Tout_model;

res = ((Ratio_model-Ratio_data)./Ratio_model).^2;

% >>>> upd. 21.11.2023: make the output human-readable: >>>> %
% Z = sqrt(sum(res))/length(res)
% X
Z = sqrt(sum(res))/length(res);
X;
fprintf("Z = %0.3e :\t",Z);
fprintf("X = ");
fprintf("%0.3e\t",X);
fprintf(".\n");
% <<<< 21.11.2023 <<<< %

figure(figNum)
clf
semilogx(tdelay, -real(Ts)./imag(Ts),'g')
hold on
semilogx(tdelay,Ratio_data,'ob',tdelay,Ratio_model,'k'); 
pause(0.01)