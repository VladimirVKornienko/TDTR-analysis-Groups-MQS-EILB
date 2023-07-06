% This program is the main program for SIMULATION of TDTR using the "Bidirectional" model.
% This program is partly based on TDTR_MAIN_V4.m published by Joseph P. Feser under http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html (12/Sept/2012).
% This program can handle varying pump size over delay time thanks to Greg Hohensee.
% This program allows you to simulate and plot TDTR signals (SIM), compare the results with measurement data (load_Data), and compute and plot sensitivity coefficients 
tic 
clear all

%-------------------------------BEGIN CODE---------------------------------

%% USER INPUT
    [SysParam] = Parameter_Data_Example();%Parameter_Example(); %load parameters (matlab function, parameters are assigned in next section below)
     load_Data = 0; if load_Data, load('struct_test_data.mat'), tnorm = 200; end   %load data (.mat), choose time value for normalization of Vin (ps)
  save_results = 0; if save_results, savefile = 'AN160219c_TDTR_Bidirectional_SENS'; end   % filename for saving parameters and results
           SIM = 1;  % do simulation
          SENS = 1;  % compute sensitivities                        
           Col = 'k';  % define color of model curve in results plot for simlation
      ClearFig = 1;  % clear results figure       
           psc = 0; if psc, frac = 0.19; tmax = 3.6e-9; end %frac: fractional change in pump spot size over delay time; tmax: maximum time delay
        nnodes = 35; %this is the number of nodes used for numerical integration. (DEFAULT=35)  
                    %35 nodes gives more than enough accuracy, even for extreme cases such as
                    %graphite with a diffraction limited spot size...
                    %(more than 5 digits of Vin,Vout precision), but if you want speed at the
                    %expense of some precision, this gives you the option. CHANGE WITH CARE!

%% LOAD PARAMETERS
% ROW vectors starting with top layer
Lambda = SysParam.Lambda; % Thermal conductivities (W m^-1 K^-1)
C = SysParam.C;  % Volumetric heat capacities (J m^-3 K^-1)
h = SysParam.h;  % Thicknesses (m)  
eta = SysParam.eta;   % Anisotropy parameter eta=kx/ky;

X_heat = SysParam.X_heat;  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
X_temp = SysParam.X_temp;  % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
AbsProf = SysParam.AbsProf;  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized

f = SysParam.f;  % Laser modulation frequence(Hz)
r_pump = SysParam.r_pump;  % Pump 1/e^2 radius (m)
r_probe = SysParam.r_pump;   % Probe 1/e^2 radius (m)
tau_rep = SysParam.tau_rep;  % Laser repetition period (s)
P_pump = SysParam.P_pump;  % absorbed pump power (transmission of objective X absorbance X pump power) 
P_probe = SysParam.P_probe;  % absorbed pump power (transmission of objective X absorbance X pump power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
tdelay_model = SysParam.tdelay_model; % time delays for model curve (s)

if psc
    r_pump_model = r_pump*(1 + frac*tdelay_model/tmax); 
else
    r_pump_model = r_pump;
end

%% IMPORT DATA
if load_Data
    tdelay_raw = Data.tdelay; % delay time (ps)
    Vin_raw = Data.Vin;  % in-phase TDTR signal (microvolts)
    Vout_raw = Data.Vout;  % out-of-phase TDTR signal (microvolts)
    Ratio_raw = Data.Ratio;  % -Vin./Vout
    Vdet_raw = Data.Vdet;  % detector voltage (mV)
end  

%% COMPUTE Tin, Tout, Ratio AND SAVE
if SIM
    [Ts,~] = TDTR_Bidirectional_SUB_B(tdelay_model,tau_rep,f,Lambda,C,h,eta,r_pump_model,r_probe,P_pump,nnodes,X_heat,X_temp);
    % Ts is a complex MATRIX with dimensions [length(tdelay_model),length(AbsProf)];
    % The ith COLUMN of Ts is the temperature response calculated at depth X_temp for heating of the sample at an infinitesimal surface at depth X_heat(i);  
    % Assuming the problem is linear, Ts is weighted with AbsProf to consider the finite heating depth described by AbsProf;
    Tin_model = real(Ts)*AbsProf./(ones(size(AbsProf))'*AbsProf);  
    Tout_model = imag(Ts)*AbsProf./(ones(size(AbsProf))'*AbsProf);
    Ratio_model = -Tin_model./Tout_model;
    if save_results
        save(strcat(pwd,'/',savefile,'.mat'));
    end
end

%% COMPUTE SENSITIVITIES
if SENS
    [S_C,S_L,S_h,S_eta,S_r_pump,S_CX,S_LX,S_hX,S_etaX,S_r_pumpX]  = TDTR_Bidirectional_SUB_SENS(tdelay_model,tau_rep,f,Lambda,C,h,eta,r_pump,r_probe,P_pump,nnodes,X_heat,X_temp,AbsProf);
    if save_results
        save(strcat(pwd,'/',savefile,'.mat'));
    end  
end

%% PLOT RESULTS
defsize = [0 0 30 15]; % set the size of the figure (left bottom width height)    
Fig = figure(111);
if ClearFig
    clf
else 
    hold on
end
set(Fig,'units','centimeters')
set(Fig,'Position',defsize)

if SIM
    subplot(1,2,1) 
    loglog(tdelay_model*1e12,Ratio_model,Col,'LineWidth',2)
    hold on
    if load_Data
        loglog(tdelay_raw,Ratio_raw,'ro','MarkerSize',4)
    end
    xlabel('time delay (s)','FontSize',18)
    ylabel('-Vin/Vout','FontSize',18)
    set(gca,'FontSize',18)
    xlim([1 5000])
    ylim([0.1 10*ceil(log10(max(Ratio_model)))])
 
    subplot(1,2,2)
    if load_Data
        Vin_norm = interp1(tdelay_raw,Vin_raw,tnorm);
        Tin_model_norm = interp1(tdelay_model*1e12,Tin_model,tnorm);
        Tin_data = Vin_raw/Vin_norm*Tin_model_norm;
        loglog(tdelay_raw,Tin_data,'ro','MarkerSize',4)
        hold on
        Vout_norm = interp1(tdelay_raw,Vout_raw,tnorm);
        Tout_model_norm = interp1(tdelay_model*1e12,Tout_model,tnorm);
        Tout_data = Vout_raw/Vout_norm*Tout_model_norm;
        loglog(tdelay_raw,-Tout_data,'ro','MarkerSize',4)        
    end
    loglog(tdelay_model*1e12,Tin_model,Col,'LineWidth',2)
    hold on
    loglog(tdelay_model*1e12,-Tout_model,Col,'LineWidth',2)
    xlabel('time delay (s)','FontSize',18)
    ylabel('\Delta T','FontSize',18)
    set(gca,'FontSize',18)
    xlim([1 5000])
    ylim([0.1 10*ceil(log10(max(Tin_model)))])
end
if SENS
    ColMat = [0 0 0; 1 0 0; 0 0 1; 0 0.5 0; 1 1 0; 0.5 0.5 0.5; 1 0.64 0];
    subplot(1,2,1) %plots subplot in top left corner of 2x2 grid; use 'position' when adding y-axis on the right side, else it does not work well
    for n = 1:length(Lambda)
        semilogx(tdelay_model*1e12,S_C(:,n),'o','Color',ColMat(n,:),'LineWidth',1,'MarkerSize',4)
        hold on
        plot(tdelay_model*1e12,S_L(:,n),'+','Color',ColMat(n,:),'LineWidth',1,'MarkerSize',4)
        plot(tdelay_model*1e12,S_h(:,n),'x','Color',ColMat(n,:),'LineWidth',1,'MarkerSize',4)
    end
    plot(tdelay_model*1e12,S_r_pump,'k--','LineWidth',1)
    if length(Lambda) == 2
        legend('S_{C(1)}','S_{L(1)}','S_{h(1)}','S_{C(2)}','S_{L(2)}','S_{h(2)}','S_r_pump','Location','NorthEastOutside');
    elseif length(Lambda) == 3
        legend('S_{C(1)}','S_{L(1)}','S_{h(1)}','S_{C(2)}','S_{L(2)}','S_{h(2)}','S_{C(3)}','S_{L(3)}','S_{h(3)}','S_r_pump','Location','NorthEastOutside');
    elseif length(Lambda) == 4
        legend('S_{C(1)}','S_{L(1)}','S_{h(1)}','S_{C(2)}','S_{L(2)}','S_{h(2)}','S_{C(3)}','S_{L(3)}','S_{h(3)}','S_{C(4)}','S_{L(4)}','S_{h(4)}','S_{rpumpprobe}','Location','NorthEastOutside');
    elseif length(Lambda) == 5
        legend('S_{C(1)}','S_{L(1)}','S_{h(1)}','S_{C(2)}','S_{L(2)}','S_{h(2)}','S_{C(3)}','S_{L(3)}','S_{h(4)}','S_{C(4)}','S_{L(4)}','S_{h(4)}','S_{C(5)}','S_{L(5)}','S_{h(5)}','S_r_pump','Location','NorthEastOutside');
    end
    set(gca,'TickLength',[0.02 0.025])
    set(gca,'LineWidth',1)
    set(gca,'FontSize',14,'FontName','Arial')
    xlim([10 5000])
    set(gca,'XTick',[10^1 10^2 10^3]);
    set(gca,'XTicklabel',[10^1 10^2 10^3]);
    % ylim([0 350])
    % set(gca,'YTick',[0.1 1 10 100])
    % set(gca,'YTicklabel',[0.1 1 10 100]);
    xlabel('\it{t}\rm (ps)','FontSize',14,'FontName','Arial')
    ylabel('Sensitivity coefficient for Ratio','FontSize',14,'FontName','Arial')
    
    subplot(1,2,2) %plots subplot in top left corner of 2x2 grid; use 'position' when adding y-axis on the right side, else it does not work well
    for n = 1:length(Lambda)
        semilogx(tdelay_model*1e12,S_CX(:,n),'o','Color',ColMat(n,:),'LineWidth',1,'MarkerSize',4)
        hold on
        plot(tdelay_model*1e12,S_LX(:,n),'+','Color',ColMat(n,:),'LineWidth',1,'MarkerSize',4)
        plot(tdelay_model*1e12,S_hX(:,n),'x','Color',ColMat(n,:),'LineWidth',1,'MarkerSize',4)
    end
    plot(tdelay_model*1e12,S_r_pumpX,'k--','LineWidth',1)
    if length(Lambda) == 2
        legend('SX_{C(1)}','SX_{L(1)}','SX_{h(1)}','SX_{C(2)}','SX_{L(2)}','SX_{h(2)}','SX_r_pump','Location','NorthEastOutside');
    elseif length(Lambda) == 3
        legend('SX_{C(1)}','SX_{L(1)}','SX_{h(1)}','SX_{C(2)}','SX_{L(2)}','SX_{h(2)}','SX_{C(3)}','SX_{L(3)}','SX_{h(3)}','SX_r_pump','Location','NorthEastOutside');
    elseif length(Lambda) == 4
        legend('SX_{C(1)}','SX_{L(1)}','SX_{h(1)}','SX_{C(2)}','SX_{L(2)}','SX_{h(2)}','SX_{C(3)}','SX_{L(3)}','SX_{h(3)}','SX_{C(4)}','SX_{L(4)}','SX_{h(4)}','SX_{rpumpprobe}','Location','NorthEastOutside');
    elseif length(Lambda) == 5
        legend('SX_{C(1)}','SX_{L(1)}','SX_{h(1)}','SX_{C(2)}','SX_{L(2)}','SX_{h(2)}','SX_{C(3)}','SX_{L(3)}','SX_{h(4)}','SX_{C(4)}','SX_{L(4)}','SX_{h(4)}','SX_{C(5)}','SX_{L(5)}','SX_{h(5)}','SX_r_pump','Location','NorthEastOutside');
    end
    set(gca,'TickLength',[0.02 0.025])
    set(gca,'LineWidth',1)
    set(gca,'FontSize',14,'FontName','Arial')
    xlim([10 5000])
    set(gca,'XTick',[10^1 10^2 10^3]);
    set(gca,'XTicklabel',[10^1 10^2 10^3]);
    % ylim([0 350])
    % set(gca,'YTick',[0.1 1 10 100])
    % set(gca,'YTicklabel',[0.1 1 10 100]);
    xlabel('\it{t}\rm (ps)','FontSize',14,'FontName','Arial')
    ylabel('Sensitivity coefficient for Vin','FontSize',14,'FontName','Arial')
end

toc
