%TDTR_matrix_to_struct
clear;
filename = ['Quartz_N4_June2.mat'];
%filename = 'f07_SampleN4_phase0_250322_230743.mat';
%filename = 'f06_SampleN3_phase0_250322_224312.mat';
%filename = 'struct_test_data2_SHIFTED.mat'
%filename = '25_f05_SampleN2_phase0_250322_220940.mat';
%filename = 'sio2_test_SHIFTED.mat';
data = load(filename);

if isstruct(data)
    structData = struct();
    structData.stagePosition = data.StagePosition;%.Data(:,1);
    structData.tdelay = data.tdelay; %.Data(:,2); % in ps
    structData.Vin = data.Vin; %.Data(:,3); %in uV
    structData.Vout = data.Vout; %.Data(:,4); % in uV
    structData.Ratio = data.ratio; %.Data(:,5);
    structData.Vdet = data.Vdet;
    
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
ylabel('V_{in}')
xlabel('t_{delay}')
subplot(3,1,2)
plot(Data.tdelay,Data.Vout)
ylabel('V_{out}')
xlabel('t_{delay}')
subplot(3,1,3)
plot(Data.tdelay,Data.Ratio)
ylabel('Ratio')
xlabel('t_{delay}')

prefix = strsplit(filename,'.');
new_file = [prefix{1} '_prep.mat'];
save(new_file,"Data")


%%Automatically run AutoSetPhase() to fix phase shift emerging from
%%instrumentation

AutoSetPhase(0,new_file,[-20 80])