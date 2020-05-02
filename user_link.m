clc;clear;close all; 
set(0,'DefaultTextInterpreter','default');
set(0, 'DefaultLegendInterpreter', 'default');


% EIRP needed
f = 4.01e8; %[Hz] %UHF
c = 3e8; %[m/s]
lambda = c / f;

R_mars = 3389.5e3; %[m]


% power_trans = [1, 5, 10, 20, 30];

min_power_recieved = -140; %dBW
gain_reciever = 12; %dBi
design_margin = 3; %dB

path_loss = @(x) 20 * log10(4 * pi * x / lambda);

SMA = linspace(100e3, 20000e3, 1000) + R_mars; %[m]


losses = path_loss(SMA - R_mars);

EIRP = min_power_recieved - gain_reciever + losses + design_margin;

figure
plot(SMA/1000, EIRP, 'LineWidth', 2)
xlabel('SMA [km]')
ylabel('EIRP [dBW]')
title('EIRP needed for minimum power')
hold on
xline(R_mars/1000);
xline(11500);
grid on

% Different Gains
gains = [0, 2, 5, 10, 15];
n = length(gains);

figure
hold on
for kk = 1 : n
    P_db = EIRP - gains(kk);
    P = 10.^(P_db./10);
    plot(SMA/1000, P, 'LineWidth', 2)
end
xlabel('SMA [km]')
ylabel('P [W]')
title('Power needed with different gains @ UHF')
hold on
xline(R_mars/1000);
xline(11500);
grid on
legend('G = 0', 'G = 2', 'G = 5', 'G = 10', 'G = 15')

% Check CRN
k = 1.3806e-23;
T_bb = 20;
T_amp = 30;

CNR_min = linspace(3, 25, 100);

figure
hold on

    P_noise = min_power_recieved - CNR_min;
    bandwidth = 10.^(P_noise./ 10) ./ (k*(T_bb + T_amp));
    plot(CNR_min, bandwidth/1e6, 'LineWidth', 2)

xlabel('Minimum CNR [dB]')
ylabel('Bandwidth [MHz]')
title('Maximum bandwidth to achieve a given CNR @UHF')
hold on
yline(4);
yline(1);
grid on









