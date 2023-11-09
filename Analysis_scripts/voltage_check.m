% First section creates 3 tiles of plots to a figure: Vin, Vout, and ratio
% Also some testing for averaging voltage values if avg_Vout_flag is true

% Second section plots 4 different avagings of ratio in a figure:
% 1 - no averaging
% 2 - moving average of Vout (currently averages 5 data points)
% 3 - moving average of Vin and Vout
% 4 - moving average of ratio

% Plot figures of the voltages
close all
clear all
%filename1 = 'dataForFits\Si_2MHz_180423_130021_FIN_NoNaNs.mat';
%filename1 = 'dataForFits\Si_050nmSiO2_2MHz_180423_161212_FIN_NoNaNs.mat';


%load(filename1)



%subplot(3,1,1)
%hold on
%plot(Data.tdelay*1e-12,Data.Vin,'b')


%subplot(3,1,2)
%hold on
%plot(Data.tdelay*1e-12,Data.Vout,'b')


%subplot(3,1,3)
%hold on
%semilogx(Data.tdelay*1e-12,Data.Ratio,'b')

%ylabel('ratio')


clear all

avg_Vout_flag = false; % Set Vout values after certain tdelay to an average in between tdelay_avg_min and tdelay_avg_max
% an average Vout is calculated from values in between tdelay_avg_min and
% tdelay_avg_max. Values after tdelay_avg_max are set to that constant
% value
tdelay_avg_min = 3; % picoseconds
tdelay_avg_max = 2000; % picoseconds

projectPath = 'C:\Users\immonel1\Downloads\Analysis_scripts';
dataFolder = '\2023_07_12 CuI new\f03_preprocessed_data';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a01_Si_Al80nm_point1_sequence1_150623_131327_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a02_Si_second_run_different_sequence_corrected_150623_141649_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a11_Si_SiO2_50nm_NDto800nmAdded_corrected_150623_200157_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a03_CuI_0deg_GLAD_150623_153912_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a04_CuI_point2_redone_150623_165645_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a05_CuI_70degr_GLAD_point1_150623_173101_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a06_CuI_70degr_GLAD_point2_corrected_150623_175319_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a07_CuI_75degr_glad_point1_150623_182217_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a08_CuI_75degr_glad_point2_corrected_150623_184625_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a09_CuI_80degr_glad_point1_150623_191225_FIN.mat';
%filename2 = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a10_CuI_80degr_glad_point2_150623_193523_FIN.mat';
%filename2 = '\f01_Si_point1_120723_153555_FIN.mat';
%filename2 = '\f02_Si_point2_120723_155842_FIN.mat';
%filename2 = '\f03_CuI_0deg_point1_120723_163922_FIN.mat';
%filename2 = '\f04_CuI_0deg_point2_120723_172648_FIN.mat';
%filename2 = '\f05_CuI_80deg_1min_point1_120723_175747_FIN.mat';
%filename2 = '\f06_CuI_80deg_1min_point2_120723_182727_FIN.mat';
%filename2 = '\f07_CuI_80deg_2min_point1_120723_185415_FIN.mat';
filename2 = '\f08_CuI_80deg_2min_point2_120723_192033_FIN.mat';
filename2 = [projectPath dataFolder filename2];

plotFontSize = 10;

load(filename2)

i = find(Data.tdelay > 0); % First positive tdelay index

tiledlayout(3,1, 'Padding', 'none', 'TileSpacing', 'compact');

%subplot(3,1,1)
nexttile
semilogx(Data.tdelay(i)*1e-12,Data.Vin(i),'b')
xlabel('\it{t}_{d} (s)')
ylabel('{\it V}_{in} (mV)')

if avg_Vout_flag
    %Vout_raw = Data.Vout;
    %j = find(Data.tdelay > tdelay_avg_min);
    %k = find(Data.tdelay > tdelay_avg_max);


    %Vout_raw(k) = mean(Data.Vout(j:k)); % i:j gives the correct indices, not sure why
    Vout_raw = average_voltage(Data.Vout,Data.tdelay,tdelay_avg_min,tdelay_avg_max);
    %4e-11
end
%fontsize(plotFontSize, 'points')
% Linear fit for Vout
%j = find(Data.tdelay > 40);
%p = polyfit(Data.tdelay(j),Data.Vout(j),1);
%V_out_fit = polyval(p,Data.tdelay(j));

%subplot(3,1,2)
nexttile

if avg_Vout_flag
    h = gobjects(2, 1);
    h(1) = semilogx(Data.tdelay(i)*1e-12,Vout_raw(i),'r','DisplayName','avg');
    hold on
    h(2) = semilogx(Data.tdelay(i)*1e-12,Data.Vout(i),'b','DisplayName','raw');
    legend(h([2 1]),'Location','northwest')
else
    semilogx(Data.tdelay(i)*1e-12,Data.Vout(i),'b','DisplayName','raw');
    
    %Linear fit for Vout
    %hold on
    %semilogx(Data.tdelay(j)*1e-12,V_out_fit,'r')
end
%fontsize(plotFontSize, 'points')
xlabel('\it{t}_{d} (s)')
ylabel('{\it V}_{out} (mV)')

%subplot(3,1,3)
nexttile


if avg_Vout_flag
    h = gobjects(2,1);
    h(1) = semilogx(Data.tdelay(i)*1e-12,-Data.Vin(i)./Vout_raw(i),'r','DisplayName','avg');
    hold on
    h(2) = semilogx(Data.tdelay(i)*1e-12,Data.Ratio(i),'b','DisplayName','raw');
    legend(h([2 1]))
else
    semilogx(Data.tdelay(i)*1e-12,Data.Ratio(i),'b','DisplayName','raw')
end
%xlim('padded')
xlabel('\it{t}_{d} (s)')
ylabel('-{\it V}_{in}/{\it V}_{out} (a.u.)')
%fontsize(plotFontSize, 'points')
if avg_Vout_flag
    tdelay_min = 0e-12; % before Apr19: 30e-12; 
    tdelay_max = 5000e-12;
    movMeanWindow = 5;
    
    tdelay_raw = Data.tdelay*1e-12; % delay time (s)
    Vin_raw = Data.Vin;  % in-phase TDTR signal (microvolts)
    %Vout_raw = Data.Vout;  % out-of-phase TDTR signal (microvolts)
    Ratio_raw = -Vin_raw./Vout_raw;  % -Vin./Vout
    Vdet_raw = Data.Vdet;  % detector voltage (mV)
    
    [~,Vin_data] = extract_interior_V4(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
    [~,Vout_data] = extract_interior_V4(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);
    [tdelay_data,Ratio_data] = extract_interior_V4(tdelay_raw,Ratio_raw,tdelay_min,tdelay_max);
    
    
    Vout_data_avg = movmean(Vout_data,movMeanWindow);
    Vin_data_avg = movmean(Vin_data,movMeanWindow);
    % average Vout
    Ratio_data_Vout = -Vin_data./Vout_data_avg;
    % average Vout and Vin
    Ratio_data_VoutVin = -Vin_data_avg./Vout_data_avg;
    % average ratio
    Ratio_data_Ratio = movmean(Ratio_data,movMeanWindow);
    % raw data as is
    
    markerSize = 2;
    figure(2)
    semilogx(tdelay_data,Ratio_data,'b-o','MarkerSize',12, 'DisplayName', 'raw');hold on
    semilogx(tdelay_data,Ratio_data_Vout,'r-d','MarkerSize',10, 'DisplayName', 'Vout')
    %hold on
    semilogx(tdelay_data,Ratio_data_VoutVin,'k-x','MarkerSize',8, 'DisplayName', 'Vout Vin')
    semilogx(tdelay_data,Ratio_data_Ratio,'g-s','MarkerSize',2, 'DisplayName', 'ratio')
    legend('Location','northeast')
    xlabel('\it{t}_{d} (s)')
    ylabel('-{\it V}_{in}/{\it V}_{out} (a.u.)')
end

%%
clear all
close all

projectPath = 'C:\Users\immonel1\Downloads\Analysis_scripts';
dataFolder = '\2023_07_12 CuI new\f03_preprocessed_data';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a01_Si_Al80nm_point1_sequence1_150623_131327_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a02_Si_second_run_different_sequence_corrected_150623_141649_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a11_Si_SiO2_50nm_NDto800nmAdded_corrected_150623_200157_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a03_CuI_0deg_GLAD_150623_153912_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a04_CuI_point2_redone_150623_165645_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a05_CuI_70degr_GLAD_point1_150623_173101_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a06_CuI_70degr_GLAD_point2_corrected_150623_175319_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a07_CuI_75degr_glad_point1_150623_182217_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a08_CuI_75degr_glad_point2_corrected_150623_184625_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a09_CuI_80degr_glad_point1_150623_191225_FIN.mat';
%filename = '2023_06_15 CuI_data\f03 (USE THIS ONE FOR FITTING) phase_corrected data\a10_CuI_80degr_glad_point2_150623_193523_FIN.mat';
%filename = 'Tarmo_data\2023_07_07\Si_2MHz_070723_162924.mat';
%filename = '\f01_Si_point1_120723_153555_FIN.mat';
%filename = '\f02_Si_point2_120723_155842_FIN.mat';
%filename = '\f03_CuI_0deg_point1_120723_163922_FIN.mat';
%filename = '\f04_CuI_0deg_point2_120723_172648_FIN.mat';
%filename = '\f05_CuI_80deg_1min_point1_120723_175747_FIN.mat';
%filename = '\f06_CuI_80deg_1min_point2_120723_182727_FIN.mat';
%filename = '\f07_CuI_80deg_2min_point1_120723_185415_FIN.mat';
filename = '\f08_CuI_80deg_2min_point2_120723_192033_FIN.mat';
filename = [projectPath dataFolder filename];
load(filename)

tdelay_min = 0e-12; % before Apr19: 30e-12; 
tdelay_max = 5000e-12;
movMeanWindow = 5;

tdelay_raw = Data.tdelay*1e-12; % delay time (s)
Vin_raw = Data.Vin;  % in-phase TDTR signal (microvolts)
Vout_raw = Data.Vout;  % out-of-phase TDTR signal (microvolts)
Ratio_raw = Data.Ratio;  % -Vin./Vout
Vdet_raw = Data.Vdet;  % detector voltage (mV)

% These lines create a linear fit for Vout after set time
%j = find(Data.tdelay > 40);
%p = polyfit(Data.tdelay(j),Data.Vout(j),1);
%V_out_fit = polyval(p,Data.tdelay(j));
%test1=Vout_raw;
%Vout_raw(j) = V_out_fit;
%test2=Vout_raw;
%Ratio_raw = -Vin_raw./Vout_raw;

[~,Vin_data] = extract_interior_V4(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
[~,Vout_data] = extract_interior_V4(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);
[tdelay_data,Ratio_data] = extract_interior_V4(tdelay_raw,Ratio_raw,tdelay_min,tdelay_max);


Vout_data_avg = movmean(Vout_data,movMeanWindow);
Vin_data_avg = movmean(Vin_data,movMeanWindow);
% average Vout
Ratio_data_Vout = -Vin_data./Vout_data_avg;
% average Vout and Vin
Ratio_data_VoutVin = -Vin_data_avg./Vout_data_avg;
% average ratio
Ratio_data_Ratio = movmean(Ratio_data,movMeanWindow);
% raw data as is

linear = 0; % 1 for a linear plot, 0 for semilogx
markerSize = 2;
figure(1)
if linear
    plot(tdelay_data,Ratio_data,'b-o','MarkerSize',markerSize, 'DisplayName', 'raw');hold on
    plot(tdelay_data,Ratio_data_Vout,'r-o','MarkerSize',markerSize, 'DisplayName', 'Vout')
    %hold on
    plot(tdelay_data,Ratio_data_VoutVin,'k-o','MarkerSize',markerSize, 'DisplayName', 'Vout Vin')
    plot(tdelay_data,Ratio_data_Ratio,'g-o','MarkerSize',markerSize, 'DisplayName', 'ratio')
else
    semilogx(tdelay_data,Ratio_data,'b-o','MarkerSize',markerSize, 'DisplayName', 'raw');hold on
    semilogx(tdelay_data,Ratio_data_Vout,'r-o','MarkerSize',markerSize, 'DisplayName', 'Vout')
    %hold on
    semilogx(tdelay_data,Ratio_data_VoutVin,'k-o','MarkerSize',markerSize, 'DisplayName', 'Vout Vin')
    semilogx(tdelay_data,Ratio_data_Ratio,'g-o','MarkerSize',markerSize, 'DisplayName', 'ratio')
    %hold on
    %plot([1.20833 1.20833]*1e-11, [0.4 0.8], 'k', 'HandleVisibility','off')
    %plot([3.70833 3.70833]*1e-11, [0.4 0.8], 'k', 'HandleVisibility','off')
    %plot([1.05 1.05]*1e-10, [0.4 0.8], 'k', 'HandleVisibility','off')
%xlim([1e-12 1e-9])
end
legend('Location','northeast')
xlabel('\it{t}_{d} (s)')
ylabel('-{\it V}_{in}/{\it V}_{out} (a.u.)')
fontsize(13,'points')

if false
    figure(2)
    semilogx(tdelay_data,Vout_data,'b-o','MarkerSize',markerSize,'DisplayName','raw');hold on
    semilogx(tdelay_data,Vout_data_avg,'r-o','MarkerSize',markerSize,'DisplayName','avg')
    legend('Location', 'northeast')
    xlabel('\it{t}_{d} (s)')
    ylabel('{\it V}_{out} (mV)')
    fontsize(13,'points')
end