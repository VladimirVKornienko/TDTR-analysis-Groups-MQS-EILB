%Pick positions of acoustic echoes in GUI and calculate the transducer film
%thickness


clear;
filename = uigetfile('*.mat');
data = load(filename);

Data = data.Data;

figure;
% subplot(3,1,1)
% plot(Data.tdelay,Data.Vin)
% ylabel('Vin')
% xlabel('t_{delay} (ps)')
% subplot(3,1,2)
% plot(Data.tdelay,Data.Vout)
% ylabel('Vout')
% xlabel('t_{delay} (ps)')
% subplot(3,1,3)
plot(Data.tdelay,Data.Ratio)
ylabel('Ratio')
xlabel('t_{delay} (ps)')
xlim([-20 150])
hold on;


%prompt = 'Enter the number of echoes to be selected:'
%n_e = input(prompt);
%[t_e, mag] = ginput(n_e);

[t_e, mag] = ginput();
title('Select echoes and press Enter')

n_e = length(t_e);

scatter(t_e,mag,'+','k')

t_echo = diff(t_e);
avg_t_echo = sum(t_echo)/(n_e-1);

v_s  = 6260; % m/s, for Al
d_transducer = 0.5*avg_t_echo*1e-12*v_s; %convert t_echo from ps to s

d_transducer_nm = d_transducer*1e9;

disp(['Transducer thickness estimated from echo: ',num2str(d_transducer_nm),' nm'])
