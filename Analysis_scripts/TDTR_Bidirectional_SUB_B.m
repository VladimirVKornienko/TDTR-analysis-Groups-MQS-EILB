% This function is part of the "Bidirectional" TDTR analysis.
% This function computes the integral in Eq. (9), the summation in Eqs. (19) and (20), and the phase shift in Eq. (21) in Ref. 2004_Cahill_RSI.
% This function is based on TDTR_REFL_V4.m published by Joseph P. Feser under http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html (12-September-2012).
% This function can handle varying pump size over delay time thanks to Greg Hohensee.

% Input parameters:
%   tdelay   - time delay (ps)
%   rau_rep  - repetition rate of laser (s)
%   f        - modulation frequency (Hz)
%   Lambda   - vector of thermal conductivities, Lambda(1) = top surface,(W/m/K)
%   C        - vector of volumetric heat capacity (J/m^3/K)
%   h        - thicknesses of each layer (layer N will NOT be used, semiinfinite)
%   r_pump   - pump 1/e^2 radius (m)
%   r_probe  - probe 1/e^2 radius (m)
%   P_pump   - absorbed pump laser power (W) ...doesn't effect ratio
%   X_heat   - temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
%   X_temp   - depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)

%---------------------------- BEGIN CODE ----------------------------------

function [T,Ratio] = TDTR_Bidirectional_SUB_B(tdelay,tau_rep,f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,X_heat,X_temp)

ii=sqrt(-1);
fmax=20/min(abs(tdelay)); %maximum frequency considered (see RSI paper) default:10/min(abs(tdelay)) for accuracy at short delay times, increase to 20
kmax = 1/sqrt(r_pump(end)^2+r_probe^2)*2; %cutoff wavevector for integration
%use Legendre-Gauss Integration
%computes weights and node locations...
if nargin<12
    nnodes = 35;
end
[kvect,weights] = lgwt_V4(nnodes,0,kmax);

M=20*ceil(tau_rep/min(abs(tdelay))); %Highest Fourier component considered
mvect=-M:M; %Range of Fourier components to consider (Vector)
fudgep=exp(-pi*((mvect/tau_rep+f)/fmax).^2);%artificial decay (see RSI paper)
fudgem=exp(-pi*((mvect/tau_rep-f)/fmax).^2);

T = zeros(length(tdelay),length(X_heat));
Ratio = T;
for i = 1:length(X_heat)
    Ip = TDTR_Bidirectional_SUB_A(kvect,mvect/tau_rep+f,Lambda,C,h,eta,r_pump,r_probe,P_pump,X_heat(i),X_temp);
    Im = TDTR_Bidirectional_SUB_A(kvect,mvect/tau_rep-f,Lambda,C,h,eta,r_pump,r_probe,P_pump,X_heat(i),X_temp);

    expterm = exp(ii*2*pi/tau_rep*(tdelay*mvect));

    [Nk,Nf,Nt] = size(Ip); % get 3D matrix dimensions
    dTp = reshape(weights' * Ip(:,:), 1, Nf, Nt);
    dTm = reshape(weights' * Im(:,:), 1, Nf, Nt); 
        
    NNt = ones(length(tdelay),1);
    if length(r_pump) > 1
        % For simplicity let's convert [1 Nf Nt] to [Nt Nf]
        % as was the original format for Retemp and Imtemp
        dTp = permute(squeeze(dTp),[2,1]);
        dTm = permute(squeeze(dTm),[2,1]);
        % NNt(Nt,1) * fudge(1,Nf) has size [Nt Nf].
        Retemp =     (dTp.*(NNt*fudgep)+dTm.*(NNt*fudgem)).*expterm;
        Imtemp = -ii*(dTp        -dTm        ).*expterm;
    else % r_pump is scalar
        % These are the original Retemp and Imtemp, for scalar r_pump
        Retemp =    (NNt*(dTp.*fudgep+dTm.*fudgem)).*expterm;
        Imtemp = -ii*(NNt*(dTp        -dTm        )).*expterm;
    end

    Resum = sum(Retemp,2); %Sum over all Fourier series components
    Imsum = sum(Imtemp,2);

    Tm = Resum+ii*Imsum; %

    %This exp term here is the phase shift fix term used to deal with "pump
    %on a stage" setup configuration. This is not the case with AALTO TDTR
    %setup and the term should thus be removed when dealing with data
    %measured with AALTO TDTR setup.
    
    T(:,i) = Tm; %.*exp(ii*2*pi*f*tdelay); %Reflectance Fluxation (Complex)

    Ratio(:,i) = -real(T(:,i))./imag(T(:,i));
end