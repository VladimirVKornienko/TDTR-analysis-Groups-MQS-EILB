%TDTR_matrix_to_struct
clear;
filename = 'f07_SampleN4_phase0_250322_230743.mat';
%filename = 'f06_SampleN3_phase0_250322_224312.mat';
%filename = 'struct_test_data2_SHIFTED.mat'
%filename = '25_f05_SampleN2_phase0_250322_220940.mat';
%filename = 'sio2_test_SHIFTED.mat';
data = load(filename);

if isstruct(data)
    structData = struct();
    structData.stagePosition = data.Data(:,1);
    structData.tdelay = data.Data(:,2); % in ps
    structData.Vin = data.Data(:,3); %in uV
    structData.Vout = data.Data(:,4); % in uV
    structData.Ratio = data.Data(:,5);
    structData.Vdet = ones(length(data.Data(:,5)),1);
    
    Data = structData;
    
elseif ismatrix(data)
    Data.stagePosition = data(:,1);
    Data.tdelay = data(:,2);
    Data.Vin = data(:,3);
    Data.Vout = data(:,4);
    Data.Ratio = data(:,5);
    Data.Vdet = ones(length(data.Data(:,5)),1);
end

%automatically shift the tdelay zero to peak position
[approx_zero,az_idx] = max(abs(diff(Data.Vin)));
Data.tdelay = (Data.stagePosition - Data.stagePosition(az_idx)).*(2*12.5e-6*1e12)./3e8; %in ps


figure;
subplot(3,1,1)
plot(Data.tdelay,Data.Vin)
ylabel('Vin')
xlabel('t_delay')
subplot(3,1,2)
plot(Data.tdelay,Data.Vout)
ylabel('Vout')
xlabel('t_delay')
subplot(3,1,3)
plot(Data.tdelay,Data.Ratio)
ylabel('Ratio')
xlabel('t_delay')

save('Quartz_N4_testdata.mat',"Data")