% This program is the main program for FITTING of the "Bidirectional" model to TDTR data.
% This program is based on TDTR_MAIN_V4.m published by Joseph P. Feser under http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html (12/Sept/2012).
% This program can handle varying pump size over delay time thanks to Greg Hohensee.

tic 
clear all

%%  Note! This fitting script accept the pre-processed data.  %%
%% The script only uses the RATIO field, if not using the smoothing (then, Vin and Vout also). %%

%% TODO later: remove smoothing into the pre-processing step! %%

%-------------------------------BEGIN CODE---------------------------------

%% USER INPUT
   
    % 20.11.2023: Use 2 pairs of parameter + data files for fitting >>>
    % 2-nd parameter file: take mod. freq. only!! %
    % Fit parameters should better coincide, but they are taken from the
    % 1st file only.
  flagTwoFiles = true;
    % <<< %

    [SysParam] = curr_Nov07_Au_on_SiO2_5MHz(); %Parameter_Example(); % load parameters (matlab function, parameters are assigned in next section below)
    if (flagTwoFiles)   % 20.11.2023 %
        [SysParam2ndFile] = curr_Nov07_Au_on_SiO2_2MHz();
    end

    
    datafile = SysParam.filename; %'Data_Example.mat';  % load data (.mat)
    if (flagTwoFiles)   % 20.11.2023 %
        f2_datafile = SysParam2ndFile.filename;
    end                 % <<<<< %
         tnorm = 200;  % choose time value for normalization of Vin (ps)
       auto_on = 1;  % 1 for automatic fitting, 0 for manual fitting
  save_results = 1; if save_results, addfilename = ''; end  % saves results using datafile string as filename; you can add further information to the filename using the addfilename string
           Col = 'k';   % Change color of curve in results plots
      ClearFig = 1;  % Clear results figure
           psc = 0; if psc, frac = 0.19; tmax = 3.6e-9; end %frac: fractional change in pump spot size over delay time; tmax: maximum time delay
        nnodes = 35; %this is the number of nodes used for numerical integration. (DEFAULT=35)  
                    %35 nodes gives more than enough accuracy, even for extreme cases such as
                    %graphite with a diffraction limited spot size...
                    %(more than 5 digits of Vin,Vout precision), but if you want speed at the
                    %expense of some precision, this gives you the option. CHANGE WITH CARE!

                    
    %                   >>>>>>>>>>>>                        %
    % OBSOLETE! Now it is done at the pre-processing stage. %
    %                                                       %
    %                                                       %
    % Updated by Vihtori: >>>>> %                           %
    %flagUseMovMean = false;   % true or false: apply the moving average transform of the size
    flagUseMovMean = 0;  % Apply the moving average transform to the given value.
    % 0: raw data (no smoothing), 1: V_out, 2: V_in and V_out, 3: Ratio.
    %                                                       %
    % movMeanWindow = 5;      % << (that size) to the RATIO data. %
    %                                                       %
    %                   <<<<<<<<<<<<                        %


%% LOAD PARAMETERS
% ROW vectors starting with top layer

% loading the 1-st file:

Lambda = SysParam.Lambda; % Thermal conductivities (W m^-1 K^-1)
C = SysParam.C;  % Volumetric heat capacities (J m^-3 K^-1)
h = SysParam.h;  % Thicknesses (m)
eta = SysParam.eta;   % Anisotropy parameter eta=kx/ky;

X_heat = SysParam.X_heat;  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated
X_temp = SysParam.X_temp;  % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
AbsProf = SysParam.AbsProf;  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized

f = SysParam.f;  % Laser modulation frequence(Hz)
r_pump = SysParam.r_pump;  % Pump 1/e^2 radius (m)
r_probe = SysParam.r_probe;   % Probe 1/e^2 radius (m)  % CORRECTED FROM ".r_pump", July 2023. %
tau_rep = SysParam.tau_rep;  % Laser repetition period (s)
P_pump = SysParam.P_pump;  % absorbed pump power (transmission of objective X absorbance X pump power)
P_probe = SysParam.P_probe;  % absorbed pump power (transmission of objective X absorbance X pump power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
tdelay_model = SysParam.tdelay_model; % time delays for model curve (s)

%  Layer indices of parameters to be adjusted in the fit
FITNLambda = SysParam.FITNLambda;
FITNC = SysParam.FITNC;
FITNh = SysParam.FITNh;
% Range of time delay range to fit (s)
tdelay_min = SysParam.tdelay_min;
tdelay_max = SysParam.tdelay_max;

if psc
    r_pump_model = r_pump*(1 + frac*tdelay_model/3.65e-9);
else
    r_pump_model = r_pump;
end


if (flagTwoFiles)  % << upd. 20.11.2023 << %

    % loading the 2-nd file:

    % JUST THE MODULATION FREQUENCY !! %
    % Also adding f2_tdelay_min and f2_tdelay_max for compatibility. %

    % f2_Lambda = SysParam2ndFile.Lambda; % Thermal conductivities (W m^-1 K^-1)
    % f2_C = SysParam2ndFile.C;  % Volumetric heat capacities (J m^-3 K^-1)
    % f2_h = SysParam2ndFile.h;  % Thicknesses (m)  
    % f2_eta = SysParam2ndFile.eta;   % Anisotropy parameter eta=kx/ky;
    % 
    % f2_X_heat = SysParam2ndFile.X_heat;  % Temperature response is calculated for each entry i of the COLUMN vector X_heat, where X_heat(i) defines the ininitesimal surface that is being heated 
    % f2_X_temp = SysParam2ndFile.X_temp;  % depth at which temperature response is calculated (SCALAR); to consider depth sensitivity of TDTR, solve for various X_temp's through the optical skin depth and weight solutions according to the depth sensitivity (assuming linear problem)
    % f2_AbsProf = SysParam2ndFile.AbsProf;  % COLUMN vector describing absorption profile through depths X_heat of the sample; does not need to be normalized
    
    f2_f = SysParam2ndFile.f;  % Laser modulation frequence(Hz)
    % f2_r_pump = SysParam2ndFile.r_pump;  % Pump 1/e^2 radius (m)
    % f2_r_probe = SysParam2ndFile.r_probe;   % Probe 1/e^2 radius (m)  % CORRECTED FROM ".r_pump", July 2023. %
    % f2_tau_rep = SysParam2ndFile.tau_rep;  % Laser repetition period (s)
    % f2_P_pump = SysParam2ndFile.P_pump;  % absorbed pump power (transmission of objective X absorbance X pump power) 
    % f2_P_probe = SysParam2ndFile.P_probe;  % absorbed pump power (transmission of objective X absorbance X pump power); assumes AF chopper is OFF!  If not, then you need to multiply the probe power by 2.
    % f2_tdelay_model = SysParam2ndFile.tdelay_model; % time delays for model curve (s)
    
    %  Layer indices of parameters to be adjusted in the fit
    % f2_FITNLambda = SysParam2ndFile.FITNLambda; 
    % f2_FITNC = SysParam2ndFile.FITNC;
    % f2_FITNh = SysParam2ndFile.FITNh;
    % Range of time delay range to fit (s)
    f2_tdelay_min = SysParam2ndFile.tdelay_min;
    f2_tdelay_max = SysParam2ndFile.tdelay_max;
    % 
    % if psc
    %     f2_r_pump_model = f2_r_pump*(1 + frac*f2_tdelay_model/3.65e-9); 
    % else
    %     f2_r_pump_model = f2_r_pump;
    % end

end

%% IMPORT DATA

% >>> 20.11.2023 >>> %
% Now, <f2_Data> holds the data from the 2nd file,
% and <Data> - from the 1st. %
if (flagTwoFiles)
    load(f2_datafile);
    f2_Data = Data;
    clear Data;
end
load(datafile);
% <<<<< %


tdelay_raw = Data.tdelay*1e-12; % delay time (s)
Vin_raw = Data.Vin;  % in-phase TDTR signal (microvolts)
Vout_raw = Data.Vout;  % out-of-phase TDTR signal (microvolts)
Ratio_raw = Data.Ratio;  % -Vin./Vout
Vdet_raw = Data.Vdet;  % detector voltage (mV)
% Extract data corresponding to time delay range to fit
[~,Vin_data] = extract_interior_V4(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
[~,Vout_data] = extract_interior_V4(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);
[tdelay_data,Ratio_data] = extract_interior_V4(tdelay_raw,Ratio_raw,tdelay_min,tdelay_max);
if psc
    r_pump_data = r_pump*(1 + frac*tdelay_data/3.65e-9); 
else
    r_pump_data = r_pump;
end


if (flagTwoFiles)   % <<< 20.11.2023 <<< %
    f2_tdelay_raw = f2_Data.tdelay*1e-12; % delay time (s)
    f2_Vin_raw = f2_Data.Vin;  % in-phase TDTR signal (microvolts)
    f2_Vout_raw = f2_Data.Vout;  % out-of-phase TDTR signal (microvolts)
    f2_Ratio_raw = f2_Data.Ratio;  % -Vin./Vout
    f2_Vdet_raw = f2_Data.Vdet;  % detector voltage (mV)
    % Extract data corresponding to time delay range to fit
    [~,f2_Vin_data] = extract_interior_V4(f2_tdelay_raw,f2_Vin_raw,f2_tdelay_min,f2_tdelay_max);
    [~,f2_Vout_data] = extract_interior_V4(f2_tdelay_raw,f2_Vout_raw,f2_tdelay_min,f2_tdelay_max);
    [f2_tdelay_data,f2_Ratio_data] = extract_interior_V4(f2_tdelay_raw,f2_Ratio_raw,f2_tdelay_min,f2_tdelay_max);
    if psc
        % NB! the same parameter used as for the file No. 1 %
        f2_r_pump_data = r_pump*(1 + frac*f2_tdelay_data/3.65e-9);
    else
        f2_r_pump_data = r_pump;
    end
end

%% Moving average for <RATIO> data, if needed:
% OBSOLETE! Now it is done at the pre-processing stage.


%% DO THE FITTING

currFigN = 10;  % handle for the figure; can be incremented or kept constant -- see below.

% define initial value(s) for fit
X0 = zeros(length(FITNLambda)+length(FITNC)+length(FITNh),1);
for i = 1:length(FITNLambda)
    X0(i) = Lambda(FITNLambda(i));   
end
for i = 1:length(FITNC)
    X0(length(FITNLambda)+i) = C(FITNC(i));  
end
for i = 1:length(FITNh)
    X0(length(FITNLambda)+length(FITNC)+i) = h(FITNh(i));
end


%if (flagTwoFiles)  % << 20.11.2023 << %
%end    
% <<< Nothing here, for all the parameters in X0 are the same for both
% files!


% automatic fitting    
if auto_on == 1      
    
    if (flagTwoFiles == false)  % << 20.11.2023 << %
        % Xsol = fminsearch(@(X) TDTR_Bidirectional_SUB_C(X,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf),X0); 
        % [Z,~] = TDTR_Bidirectional_SUB_C(Xsol,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf);
        Xsol = fminsearch(@(X) TDTR_Bidirectional_SUB_C(currFigN, X,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf),X0); 
        [Z,~] = TDTR_Bidirectional_SUB_C(currFigN, Xsol,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf);
    else
        % only change for the second call of "_SUB_C": f2_f (modulation
        % frequency).
        Xsol = fminsearch(@(X) (TDTR_Bidirectional_SUB_C(currFigN, X,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf) ...
           + TDTR_Bidirectional_SUB_C(currFigN+1, X,f2_Ratio_data,f2_tdelay_data,tau_rep,f2_f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf)) ,X0); 
        
        Z = (TDTR_Bidirectional_SUB_C(currFigN, Xsol,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf) ...
            + TDTR_Bidirectional_SUB_C(currFigN+1, Xsol,f2_Ratio_data,f2_tdelay_data,tau_rep,f2_f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf));
    end

    fprintf('Data fit completed\n')
else    

% manual fitting        
    globAnsw = 1;

    while globAnsw == 1
        if (flagTwoFiles == false)  % << 20.11.2023 << %
            %TDTR_Bidirectional_SUB_C(X0,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf)
            TDTR_Bidirectional_SUB_C(currFigN,X0,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf)
        else
            TDTR_Bidirectional_SUB_C(currFigN,X0,Ratio_data,tdelay_data,tau_rep,f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf)
            TDTR_Bidirectional_SUB_C(currFigN+1,X0,f2_Ratio_data,f2_tdelay_data,tau_rep,f2_f,Lambda,C,h,eta,r_pump_data,r_probe,P_pump,nnodes,FITNLambda,FITNC,FITNh,X_heat,X_temp,AbsProf)
        end    
        currFigN = currFigN + 2;
        globAnsw = 0;
        hold on
        answ = input('Continue with thermal conductivity (Lambda)? (Yes: "1", No: "0")');
        if answ == 1
            globAnsw = 1;
            for i = 1:length(FITNLambda)
                % num = sprintf('\t%f',FITNLambda(i));
                num = sprintf('\t%e',X0(i));
                text = strcat('Input new value for lambda',num,' : > ');
                X0(i) = input(text);
            end
        end
        
        answ = input('Continue with heat capacities? (Yes: "1", No: "0")');
        if answ == 1
            globAnsw = 1;
            for i = 1:length(FITNC)
                % num = sprintf('\t%f',FITNC(i));
                num = sprintf('\t%e',X0(i+length(FITNLambda)));
                text = strcat('Input new value for C',num,' : > ');
                X0(i+length(FITNLambda)) = input(text);
            end
        end

        answ = input('Continue with film thicknesses? (Yes: "1", No: "0")');
        if answ == 1
            globAnsw = 1;
            for i = 1:length(FITNh)
                % num = sprintf('\t%f',FITNh(i));
                num = sprintf('\t%e',X0(i+length(FITNLambda)+length(FITNC)));
                text = strcat('Input new value for h',num,' : > ');
                X0(i+length(FITNLambda)+length(FITNC)) = input(text);
            end
        end

        % PUMP POWER :: it does not affect the fit...
        % answ = input('Continue with pump power? (Yes: "1", No: "0")');
        % if answ == 1
        %    globAnsw = 1;
        %    
        %    num = sprintf('\t%e\t',P_pump);
        %    text = strcat('Input new value for Ppump',num,' : > ');
        %    P_pump = input(text);
        %    
        % end

    end
    Xsol = X0;
end

% assign fit results to adjusted parameters
for i = 1:length(FITNLambda)
    Lambda(FITNLambda(i)) = Xsol(i);   
end
for i = 1:length(FITNC)
    C(FITNC(i)) = Xsol(length(FITNLambda)+i);   
end
for i = 1:length(FITNh)
    h(FITNh(i)) = Xsol(length(FITNLambda)+length(FITNC)+i);
end

fprintf("Fit results are as follows (<Xsol> variable):\n");
fprintf("%3.5e\n",Xsol);

% temp >>> ;
return
% <<< tmp ;

%% COMPUTE Tin, Tout, Ratio AND SAVE PARAMETERS
[Ts,~] = TDTR_Bidirectional_SUB_B(tdelay_model,tau_rep,f,Lambda,C,h,eta,r_pump_model,r_probe,P_pump,nnodes,X_heat,X_temp);

Tin_model = real(Ts)*AbsProf./(ones(size(AbsProf))'*AbsProf);
Tout_model = imag(Ts)*AbsProf./(ones(size(AbsProf))'*AbsProf);
Ratio_model = -Tin_model./Tout_model;

if save_results
    save(strcat(pwd,'/',datafile(1:end-4),'_FITresults',addfilename,'.mat'));
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

subplot(1,2,1) 
loglog(tdelay_raw*1e12,Ratio_raw,'ro','MarkerSize',4)
hold on
loglog(tdelay_model*1e12,Ratio_model,Col,'LineWidth',2)
xlabel('time delay (s)','FontSize',18)
ylabel('-Vin/Vout','FontSize',18)
set(gca,'FontSize',18)
xlim([1 5000])
ylim([0.1 10*ceil(log10(max(Ratio_model)))])

subplot(1,2,2)        
loglog(tdelay_model*1e12,Tin_model,Col,'LineWidth',2)
hold on
loglog(tdelay_model*1e12,-Tout_model,Col,'LineWidth',2)
Vin_norm = interp1(tdelay_raw*1e12,Vin_raw,tnorm);
Tin_model_norm = interp1(tdelay_model*1e12,Tin_model,tnorm);
Tin_data = Vin_raw/Vin_norm*Tin_model_norm;
Vout_norm = interp1(tdelay_raw*1e12,Vout_raw,tnorm);
Tout_model_norm = interp1(tdelay_model*1e12,Tout_model,tnorm);
Tout_data = Vout_raw/Vout_norm*Tout_model_norm;
loglog(tdelay_raw*1e12,-Tout_data,'ro','MarkerSize',4)
loglog(tdelay_raw*1e12,Tin_data,'ro','MarkerSize',4)
xlabel('time delay (s)','FontSize',18)
ylabel('\Delta T','FontSize',18)
set(gca,'FontSize',18)
xlim([1 5000])
ylim([0.1 10*ceil(log10(max(Tin_model)))])

toc