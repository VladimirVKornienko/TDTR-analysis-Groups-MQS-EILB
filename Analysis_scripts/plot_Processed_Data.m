%Plot processed data (such as example or that coming from AutoSetPhase)


%TDTR_matrix_to_struct
clear;
filename = ['Al_Si@9_prep.mat'];
data = load(filename);

Data = data.Data;

figure;
subplot(3,1,1)
plot(Data.tdelay,Data.Vin)
ylabel('Vin')
xlabel('t_{delay} (ps)')
subplot(3,1,2)
plot(Data.tdelay,Data.Vout)
ylabel('Vout')
xlabel('t_{delay} (ps)')
subplot(3,1,3)
plot(Data.tdelay,Data.Ratio)
ylabel('Ratio')
xlabel('t_{delay} (ps)')


% figure;
% plot(Data.tdelay,Data.Ratio,'linewidth',2)
% ylabel('Ratio (-V_{in}/V_{out})')
% xlabel('t_{delay} (ps)')
% 
% hold on;
% 
% filename = ['Si_N3_June2_prep_SHIFTED.mat'];
% data = load(filename);
% 
% Data = data.Data;
% 
% plot(Data.tdelay,Data.Ratio,'linewidth',2)
% ylabel('Ratio (-V_{in}/V_{out})')
% xlabel('t_{delay} (ps)')
% 
% hold on;
% 
% filename = ['Si_SiO2_N2_June2_prep_SHIFTED.mat'];
% data = load(filename);
% 
% Data = data.Data;
% 
% plot(Data.tdelay,Data.Ratio,'linewidth',2)
% ylabel('Ratio (-V_{in}/V_{out})')
% xlabel('t_{delay} (ps)')
% 
% legend ('SiO2', 'Si', 'Si+1um SiO2','box','off')