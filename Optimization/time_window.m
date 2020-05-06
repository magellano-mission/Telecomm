%% Time window computation
clear, clc, close all

% 1 Day capacity
capacity = 32; %[GB]
C = capacity * 8 * 1e3; % [Mb]
C = 10 * 1e3;
datarate = linspace(5, 20, 1000);

time = C ./ datarate;

time_hours = time / 3600;
figure
plot(datarate, time_hours, 'LineWidth', 2)
hold on
yline(24)
xlabel('Datarate [Mbps]')
ylabel('Dumping time [h]')
title('Dumping time needed for 100GB @ given data-rate')
grid on