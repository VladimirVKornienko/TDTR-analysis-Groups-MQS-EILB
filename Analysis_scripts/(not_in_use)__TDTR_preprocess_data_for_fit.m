%TDTR_matrix_to_struct
clear;
flagUseAutoTimeShifting = false;
manualTimeShiftingZeroPos = -11900; % in steps.

referenceVoltage = 500; % in mV >> for the normalization by the Vdet.
flagUseVoltageVDETnormalization = true;

filename = ['a03_row2_75nm_col2_350nm_FullStep250_161222_152656.mat'];
%filename = 'Nov21_OldSi_acoustics_search_range2_try2_211122_114926.mat';

%filename = 'f07_SampleN4_phase0_250322_230743.mat';
%filename = 'f06_SampleN3_phase0_250322_224312.mat';
%filename = 'struct_test_data2_SHIFTED.mat'
%filename = '25_f05_SampleN2_phase0_250322_220940.mat';
%filename = 'sio2_test_SHIFTED.mat';
data2 = load(filename);

if isstruct(data2)
    data=data2.Data;
    structData = struct();
    structData.stagePosition = data.StagePosition;%.Data(:,1);
    structData.tdelay = data.tdelay; %.Data(:,2); % in ps
    structData.Vin = data.Vin; %.Data(:,3); %in uV
    structData.Vout = data.Vout; %.Data(:,4); % in uV
    structData.Ratio = data.Ratio; %.Data(:,5);
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


% figure;
% subplot(3,1,1)
% plot(Data.tdelay,Data.Vin)
% ylabel('V_{in}')
% xlabel('t_{delay}')
% subplot(3,1,2)
% plot(Data.tdelay,Data.Vout)
% ylabel('V_{out}')
% xlabel('t_{delay}')
% subplot(3,1,3)
% plot(Data.tdelay,Data.Ratio)
% ylabel('Ratio')
% xlabel('t_{delay}')
% title('Before shifting the zero time delay:')

%% ok, block currently... %%

%%% ignore NaN data points:
%fprintf("Total 'NaN' elements (by <Ratio>): %d\n",sum(isnan(Data.Ratio)))
%%% dataNoNans = data(~isnan(data));
%%% Data = Data(:,~isnan(Data.Ratio));
%fprintf("Length before deleting: %d\n",length(Data.Ratio));
%%% Data([isnan(Data.Ratio)])=[];
%Data([isnan(Data.Ratio)])=[];
%fprintf("Length after deleting: %d\n",length(Data.Ratio));

%automatically shift the tdelay zero to peak position
if (flagUseAutoTimeShifting)
    [approx_zero,az_idx] = max(abs(diff(Data.Vin)));
    Data.tdelay = (Data.stagePosition - Data.stagePosition(az_idx)).*(2*12.5e-6*1e12)./3e8; %in ps
    fprintf("Zero time delay set corresponding to %f steps.\n",Data.stagePosition(az_idx));
else
    Data.tdelay = (Data.stagePosition - manualTimeShiftingZeroPos).*(2*12.5e-6*1e12)./3e8; %in ps
    fprintf("Time values not shifted.\n");
end

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
% title('After shifting the zero time delay:')

prefix = strsplit(filename,'.');
new_file = [prefix{1} '_prep.mat'];
save(new_file,"Data");

%%% VKorn:
filename_short = (strsplit(filename, '.mat')); % take filename w/o extension
filename_short = string(filename_short(1));

new_name = filename_short+'_VKorn.dat'
% new_name2 = filename_short+'_VKornTEST.dat'

% save(new_name,"data")

%%% New part: normalize voltages by Vdet:
if (flagUseVoltageVDETnormalization)
    myVoltsNorm = [Data.Vin Data.Vout] ./ (Data.Vdet / referenceVoltage);
else
    myVoltsNorm = [Data.Vin Data.Vout] ./ Data.Vdet;
end

myVoltsNorm = [Data.stagePosition Data.tdelay Data.Vin Data.Vout Data.Ratio Data.Vdet myVoltsNorm(:,1) myVoltsNorm(:,2)];

%writetable(struct2table(data), new_name, 'Delimiter','\t');

myVoltsNorm = ["stagePos" "tDelayPS" "Vin" "Vout" "ratio" "Vdet" "VinNORM" "VoutNORM" ; myVoltsNorm];
% tmpqqq = struct2matrix(Data);
% tmpqqq = cell2mat(struct2cell(Data));

% myVoltsNorm = [tmpqqq myVoltsNorm];

writematrix(myVoltsNorm, new_name, 'Delimiter','\t');

%%% <<<<<

%%Automatically run AutoSetPhase() to fix phase shift emerging from
%%instrumentation

[delphase,phase,fitparam] = AutoSetPhase(0,new_file,[-20 80]);
fprintf("Phase shift applied: %5.3f",phase)