% NAVIGATION SYSTEM DOWNLINK
clear, clc, close all

% Fixing the minimum SNR (from GPS)
SNR = 15; % [dB]

% Reciever Processing Gain (to extract the C/A code)
gain_proc = 30; %[dB]

% Reciever Bandwidth
bandwidth = 4e6; %[Hz]

f = 4.01e8; %[Hz] %UHF
c = 3e8; %[m/s]
lambda = c / f;

R_mars = 3389.5e3; %[m]

path_loss = @(x) 20 * log10(4 * pi * x / lambda);
SMA = linspace(100e3, 20000e3, 1000) + R_mars; %[m]
losses = path_loss(SMA - R_mars);

% Noise 
k = 1.3806e-23;
T_bb = 250;
T_amp = 30;

power_noise = 10 * log10( k * (T_bb + T_amp) * bandwidth );

CNR = SNR - gain_proc;

power_recieved = power_noise + CNR;

gain_reciever = 0; %dBi
design_margin = 3; %dB

EIRP = power_recieved - gain_reciever + losses + design_margin;

figure
plot(SMA/1000, EIRP, 'LineWidth', 2)
xlabel('SMA [km]')
ylabel('EIRP [dBW]')
title('EIRP needed for minimum power')
hold on
xline(R_mars/1000);
xline(11500);
xlim([0, 12000])
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
xlim([0, 12000])
grid on
legend('G = 0', 'G = 2', 'G = 5', 'G = 10', 'G = 15')
