% This function is part of the "Bidirectional" TDTR analysis.
% This function computes sensitivities of the ratio and of the temperature response calculated with the bidirectional TDTR model to changes in the parameters of the model.
% This function is based on TDTR_SensPlots_V4 published by Joseph P. Feser under http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html (12-September-2012).
% This function can handle varying pump size over delay time thanks to Greg Hohensee, however, I think this will not much influence the sensitivities. To save computation time I recommend to assume constant r_pump.

% changes: 
% 1/8/2017: corrected error in the calculation of S_h(ii); prior version did not consider thickness change h_temp(ii) = h(ii)*1.01, i.e., X_temp was modified for values of X_temp larger than sum(h(1:ii))

function [S_C,S_L,S_h,S_eta,S_r_pump,S_CX,S_LX,S_hX,S_etaX,S_r_pumpX] = TDTR_Bidirectional_SUB_SENS(tdelay,tau_rep,f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,X_heat,X_temp,AbsProf)
tic
% compute deltaT and ratio from parameters assumed to be known
tic
[T0,~] = TDTR_Bidirectional_SUB_B(tdelay,tau_rep,f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,X_heat,X_temp);
toc
T0 = T0*AbsProf./(ones(size(AbsProf))'*AbsProf);
ratio0 = -real(T0)./imag(T0);

%% initialize parameters
S_C = zeros(length(tdelay),length(Lambda));
S_CX = S_C;
S_L = S_C;
S_LX = S_C;
S_h = S_C;
S_hX = S_C;
S_eta = S_C;
S_etaX = S_C;

%% compute sensitivities relative to base parameters 
for ii=1:length(Lambda)
    %-------------Specific heat-----------------
    Ctemp = C;
    Ctemp(ii) = C(ii)*1.01;
    [T_C,~] = TDTR_Bidirectional_SUB_B(tdelay,tau_rep,f,Lambda,Ctemp,h,eta,r_pump,r_probe,P_pump,nnodes,X_heat,X_temp);    
    T_C = T_C*AbsProf./(ones(size(AbsProf))'*AbsProf);
    ratio_C = -real(T_C)./imag(T_C);
    delta_C = ratio_C-ratio0;
    Num = log(ratio_C)-log(ratio0);
    Denom = log(Ctemp(ii))-log(C(ii));
    S_C(:,ii) = Num/Denom;
    NumX = log(real(T_C))-log(real(T0));
    DenomX = log(Ctemp(ii))-log(C(ii));
    S_CX(:,ii) = NumX/DenomX;
    
    %-------------Thermal Conductivity----------
    Lambdatemp = Lambda;
    Lambdatemp(ii) = Lambda(ii)*1.01;
    [T_L,~] = TDTR_Bidirectional_SUB_B(tdelay,tau_rep,f,Lambdatemp,C,h,eta,r_pump,r_probe,P_pump,nnodes,X_heat,X_temp);
    T_L = T_L*AbsProf./(ones(size(AbsProf))'*AbsProf);
    ratio_L = -real(T_L)./imag(T_L);  
    delta_L = ratio_L-ratio0;
    Num = log(ratio_L)-log(ratio0);
    Denom = log(Lambdatemp(ii))-log(Lambda(ii));
    S_L(:,ii) = Num/Denom;
    NumX = log(real(T_L))-log(real(T0));
    DenomX = log(Lambdatemp(ii))-log(Lambda(ii));
    S_LX(:,ii) = NumX/DenomX;    
    
    %-------------Layer Thickess---------------
    htemp = h;
    htemp(ii) = h(ii)*1.01;
    hsum = cumsum(htemp); 
    Vtemp = ceil(hsum/X_temp);
    VtempB = Vtemp == 1;
    I_temp = sum(VtempB)+1;  % index of layer that includes X_temp 
    if ii<I_temp  % consider thickness change to maintain the correct X_temp
        X_temptemp = X_temp + htemp(ii)-h(ii);
    else 
        X_temptemp = X_temp;
    end
    [T_h,~] = TDTR_Bidirectional_SUB_B(tdelay,tau_rep,f,Lambda,C,htemp,eta,r_pump,r_probe,P_pump,nnodes,X_heat,X_temptemp);
    T_h = T_h*AbsProf./(ones(size(AbsProf))'*AbsProf);
    ratio_h = -real(T_h)./imag(T_h);    
    delta_h = ratio_h-ratio0;
    Num = log(ratio_h)-log(ratio0);
    Denom = log(htemp(ii))-log(h(ii));
    S_h(:,ii) = Num/Denom;
    NumX = log(real(T_h))-log(real(T0));
    DenomX = log(htemp(ii))-log(h(ii));
    S_hX(:,ii) = NumX/DenomX;    
    
    %-------------Anisotropy---------------
    etatemp = eta;
    etatemp(ii) = etatemp(ii)*1.01;
    [T_eta,~] =TDTR_Bidirectional_SUB_B(tdelay,tau_rep,f,Lambda,C,h,etatemp,r_pump,r_probe,P_pump,nnodes,X_heat,X_temp);
    T_eta = T_eta*AbsProf./(ones(size(AbsProf))'*AbsProf);
    ratio_eta = -real(T_eta)./imag(T_eta);    
    delta_eta = ratio_eta-ratio0;
    Num = log(ratio_eta)-log(ratio0);
    Denom = log(etatemp(ii))-log(eta(ii));
    S_eta(:,ii) = Num/Denom;
    NumX = log(real(T_eta))-log(real(T0));
    DenomX = log(etatemp(ii))-log(eta(ii));
    S_etaX(:,ii) = NumX/DenomX;    
end
 
r_pumptemp = r_pump*1.01;
r_probetemp = r_probe*1.01;
[deltaT_r_pump,~] = TDTR_Bidirectional_SUB_B(tdelay,tau_rep,f,Lambda,C,h,eta,r_pumptemp,r_probetemp,P_pump,nnodes,X_heat,X_temp);
deltaT_r_pump = deltaT_r_pump*AbsProf./(ones(size(AbsProf))'*AbsProf);
ratio_r_pump = -real(deltaT_r_pump)./imag(deltaT_r_pump);
delta_r_pump = ratio_r_pump-ratio0;
Num = log(ratio_r_pump)-log(ratio0);
Denom = log(r_pumptemp)-log(r_pump);
S_r_pump = Num/Denom;
NumX = log(real(deltaT_r_pump))-log(real(T0));
DenomX = log(r_pumptemp)-log(r_pump);
S_r_pumpX = NumX/DenomX;