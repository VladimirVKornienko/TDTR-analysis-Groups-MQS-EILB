% This function is part of the "Bidirectional" TDTR analysis.
% This function computes the argument of the integral in Eq. (9) in Ref. 2004_Cahill_RSI.
% This function is based on TDTR_TEMP_V4.m published by Joseph P. Feser under http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html (12-September-2012).
% This function can handle varying pump size over delay time thanks to Greg Hohensee.

% Note that the temperature rise actually measured by the lock-in requires factor 2/pi for P_pump. Without this factor, the temperature rise at short time delays 
% calculated equals the per-pulse heating calculated with pump fluence averaged by probe pulse.

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
%   X_heat   - X_heat defines the ininitesimal surface that is being heated 
%   X_temp   - depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)

% Changes:
% 1/7/2017: corrected error in calculation of transfer matrix for up layers; the error was in the argument of the expterm: exp(unplus*h(n+1)) instead of the correct exp(un*h(n)); the error affects modelling for X_heat > h(1) and/or X_temp > h(1), i.e., temperature measurements in layers below top layer    
%           corrected condition on line 62; prior version was wrong for X_temp > X_heat, as it assigned the wrong layer index to the X_temp, i.e., X_temp was determined for one layer above the correct layer

%---------------------------- BEGIN CODE ----------------------------------

function [Integrand,G] = TDTR_Bidirectional_SUB_A(kvectin,freq,Lambda,C,h,eta,r_pump,r_probe,P_pump,X_heat,X_temp)

hsum = cumsum(h); 

Vheat = ceil(hsum/X_heat);
VheatB = Vheat == 1;
I_heatA = sum(VheatB)+1;  % index of heating layer
if I_heatA == 1
    X_heatL = X_heat;  
else
    X_heatL = X_heat - hsum(I_heatA-1);  % convert to layer coordinate
end
test = 0;
if X_heatL == 0
    I_heat = I_heatA;
elseif X_heatL == h(I_heatA)
    I_heat = I_heatA+1;
else 
    Lambda = [Lambda(1:I_heatA) Lambda(I_heatA:end)]; % divide heating layer to consider heating via boundary condition
    C = [C(1:I_heatA) C(I_heatA:end)];
    h = [h(1:I_heatA-1) X_heatL (h(I_heatA)-X_heatL) h(I_heatA+1:end)]; 
    eta = [eta(1:I_heatA) eta(I_heatA:end)];
    I_heat = I_heatA+1;
    test = 1;
end
    
Vtemp = ceil(hsum/X_temp);
VtempB = Vtemp == 1;
I_temp = sum(VtempB)+1;  % index of layer that includes X_temp 
if I_temp == 1
    X_tempL = X_temp;  
else
    X_tempL = X_temp - hsum(I_temp-1);  % convert to layer coordinate
end
if X_temp >= X_heat && test == 1
    I_temp = I_temp+1;  % if heating layer is split and X_temp falls in second half, then increase I_temp by 1
    X_tempL = X_tempL-X_heatL;  % adapt the layer coordinate
end    
Nfreq = length(freq);
kvect = kvectin(:)*ones(1,Nfreq);
Nlayers = length(Lambda); %# of layers
Nint = length(kvectin); %# of different spatial frequencies to calculate for

%k is a COLUMN vector (actually a matrix that changes down the rows)
%f is a ROW vector

%% apply trasfer matrix for down layers
ii = sqrt(-1);
D = Lambda./C;
omega = 2*pi*freq;
q2 = ones(Nint,1)*(ii*omega./D(Nlayers));
kvect2 = kvect.^2;
un = sqrt(4*pi^2*eta(Nlayers)*kvect2+q2);
gamman = Lambda(Nlayers)*un;
Bplus = zeros(Nint,Nfreq);  % semi-infinite boundary condition for bottom layer
Bminus = ones(Nint,Nfreq);
kterm2 = 4*pi^2*kvect2;
for n = Nlayers:-1:I_heat+1
    q2 = ones(Nint,1)*(ii*omega./D(n-1));
    unminus = sqrt(eta(n-1)*kterm2+q2);
    gammanminus = Lambda(n-1)*unminus;
    AA = gammanminus+gamman;
    BB = gammanminus-gamman;
    temp1 = AA.*Bplus+BB.*Bminus;
    temp2 = BB.*Bplus+AA.*Bminus;
    expterm = exp(unminus*h(n-1));
    Bplus = (0.5./(gammanminus.*expterm)).*temp1;
    Bminus = 0.5./(gammanminus).*expterm.*temp2;
    % These next 3 lines fix a numerical stability issue if one of the
    % layers is very thick or resistive;
    penetration_logic = logical(h(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
    Bplus(penetration_logic) = 0;
    Bminus(penetration_logic) = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamman = gammanminus;
end
gamma_Iheat = gamman;
alpha_down = Bplus;
beta_down = Bminus;

%% apply transfer matrix for up layers
q2 = ones(Nint,1)*(ii*omega./D(1));
kvect2 = kvect.^2;
un = sqrt(4*pi^2*eta(1)*kvect2+q2);
gamman = Lambda(1)*un;
Bplus = exp(-2*un*h(1)); % adiabatic boundary condition for top layer (if semi-infinite, change to 'zeros'
Bminus = ones(Nint,Nfreq);
kterm2 = 4*pi^2*kvect2;
for n = 1:1:I_heat-2
    q2 = ones(Nint,1)*(ii*omega./D(n+1));
    unplus = sqrt(eta(n+1)*kterm2+q2);
    gammanplus = Lambda(n+1)*unplus;
    AA = gammanplus+gamman;
    BB = gammanplus-gamman;
    temp1 = AA.*Bplus+BB.*Bminus;
    temp2 = BB.*Bplus+AA.*Bminus;
    expterm = exp(un*h(n));
    Bplus = (0.5./(gammanplus.*expterm)).*temp1;
    Bminus = 0.5./(gammanplus).*expterm.*temp2;
    % These next 3 lines fix a numerical stability issue if one of the
    % layers is very thick or resistive;
    penetration_logic = logical(h(n+1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
    Bplus(penetration_logic) = 0;
    Bminus(penetration_logic) = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamman = gammanplus;
    un = unplus;
end
if I_heat == 1
    gamma_Iheatminus = NaN;
    alpha_up = NaN;
    beta_up = NaN;
    fprintf('Program does not work for X_heat = 0! Use "Surfaceheating" model instead!!')
else
    q2 = ones(Nint,1)*(ii*omega./D(I_heat-1));
    un = sqrt(4*pi^2*eta(I_heat-1)*kvect2+q2);
    gamma_Iheatminus = Lambda(I_heat-1)*un;
    alpha_up = Bplus;
    beta_up = Bminus;
end

%% calculate B1 and BN
BNminus = -(alpha_up+beta_up)./(gamma_Iheat.*(alpha_down - beta_down).*(alpha_up + beta_up) + gamma_Iheatminus.*(alpha_down + beta_down).*(alpha_up - beta_up));
B1minus = (alpha_down+beta_down).*BNminus./(alpha_up+beta_up);
q2 = ones(Nint,1)*(ii*omega./D(1));
u1 = sqrt(4*pi^2*eta(1)*kvect2+q2);
B1plus = B1minus.*exp(-2*u1*h(1));

%% Get alpha and beta for temperature sensing layer
if I_temp == 1
    BTplus = B1plus;
    BTminus = B1minus;
elseif I_temp == Nlayers
    BTplus = 0;
    BTminus = BNminus;
elseif I_temp < I_heat
    q2 = ones(Nint,1)*(ii*omega./D(1));
    kvect2 = kvect.^2;
    un = sqrt(4*pi^2*eta(1)*kvect2+q2);
    gamman = Lambda(1)*un;
    Bplus = exp(-2*un*h(1)); % adiabatic boundary condition for top layer (if semi-infinite, change to 'zeros'
    Bminus = ones(Nint,Nfreq);
    kterm2 = 4*pi^2*kvect2;
    for n = 1:1:I_temp-1
        q2 = ones(Nint,1)*(ii*omega./D(n+1));
        unplus = sqrt(eta(n+1)*kterm2+q2);
        gammanplus = Lambda(n+1)*unplus;
        AA = gammanplus+gamman;
        BB = gammanplus-gamman;
        temp1 = AA.*Bplus+BB.*Bminus;
        temp2 = BB.*Bplus+AA.*Bminus;
        expterm = exp(un*h(n));
        Bplus = (0.5./(gammanplus.*expterm)).*temp1;
        Bminus = 0.5./(gammanplus).*expterm.*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic = logical(h(n+1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        Bplus(penetration_logic) = 0;
        Bminus(penetration_logic) = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gamman = gammanplus;
        un = unplus;
    end        
    BTplus = Bplus.*B1minus;
    BTminus = Bminus.*B1minus;    
elseif I_temp >= I_heat
    q2 = ones(Nint,1)*(ii*omega./D(Nlayers));
    un = sqrt(4*pi^2*eta(Nlayers)*kvect2+q2);
    gamman = Lambda(Nlayers)*un;
    Bplus = zeros(Nint,Nfreq);
    Bminus = ones(Nint,Nfreq);
    kterm2 = 4*pi^2*kvect2;
    for n = Nlayers:-1:(I_temp+1)
        q2 = ones(Nint,1)*(ii*omega./D(n-1));
        unminus = sqrt(eta(n-1)*kterm2+q2);
        gammanminus = Lambda(n-1)*unminus;
        AA = gammanminus+gamman;
        BB = gammanminus-gamman;
        temp1 = AA.*Bplus+BB.*Bminus;
        temp2 = BB.*Bplus+AA.*Bminus;
        expterm = exp(unminus*h(n-1));
        Bplus = (0.5./(gammanminus.*expterm)).*temp1;
        Bminus = 0.5./(gammanminus).*expterm.*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic = logical(h(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        Bplus(penetration_logic) = 0;
        Bminus(penetration_logic) = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gamman = gammanminus;
    end
    BTplus = Bplus.*BNminus;
    BTminus = Bminus.*BNminus;
end

q2 = ones(Nint,1)*(ii*omega./D(I_temp));
un = sqrt(4*pi^2*eta(I_temp)*kvect2+q2);
G = BTplus.*exp(un*X_tempL) + BTminus.*exp(-un*X_tempL); %The layer G(k)

[Nk,Nf] = size(kvect); Nt = length(r_pump); % define dimensions
arg1 = -pi^2*(r_pump.^2+r_probe^2)/2; % column vector C, size(Nt,1)
arg2 = kvect.^2; % matrix B, size(Nk,Nf)

expterm = exp(reshape(arg2(:) * arg1', Nk, Nf, [])); % [Nk Nf Nt] 3D matrix, Aijk = Bij*Ck
Kernal = 2*pi*P_pump*(expterm .* repmat(kvect, [1 1 Nt])); % repmat projects kvect into t-space
% if memory is an issue, user can sacrifice readability to sandwich
% Kernal and expterm into the Integrand assignment. Shouldn't be necessary,
% since Kernal and expterm vanish after this function ends.

Integrand = repmat(G, [1 1 Nt]) .* Kernal; % this reduces to G.*Kernal for Nt = 1.
  
