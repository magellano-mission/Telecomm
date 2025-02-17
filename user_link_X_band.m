% UPLINK
clc;clear;close all; 

% EIRP needed
f = 8.44e8; %[Hz] %X-Band
% f = 4.01e8; %[Hz] %UHF
c = 3e8; %[m/s]
lambda = c / f;

R_mars = 3389.5e3; %[m]

design_margin = 3; %dB

path_loss = @(x) 20 * log10(4 * pi * x / lambda);

SMA = linspace(100e3, 20000e3, 1000) + R_mars; %[m]


losses = path_loss(SMA - R_mars);

k = 1.3806e-23;
T_bb = 20 + 273;
T_amp = 30;

SNR = 10^(14/10); %[]
C = 2e6; %bps
B = C / log2(1 + SNR);

N = 1; %number of users per satellite

bandwidth = N * B;
P_noise = 10 * log10(k * (T_bb + T_amp) * bandwidth);

CNR_min = 12;

P_rec = CNR_min + P_noise;
P_recieved_watts = 10^(P_rec/10);

P_sender = 10 * log10([10 15 20 35]); %dBW
gain_sender = 5; %dBi

EIRP_sender = P_sender + gain_sender;

figure
hold on
for kk = 1 : 4
    gain_reciever = P_rec - EIRP_sender(kk) + losses + design_margin;
    plot(SMA/1000, gain_reciever, 'LineWidth', 2)
end
xlabel('SMA [km]')
ylabel('Reciever gain [dB]')
title('Reciever gain needed with CNR = 12 @ X-Band')
hold on
xline(R_mars/1000);
xline(11500);
yline(5); %60deg
yline(6); %40deg
yline(10);%20deg
grid on
legend('P = 10 W', 'P = 15 W', 'P = 20 W')







