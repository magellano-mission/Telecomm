%% Noise Power vs bandwidth frequency 100MHZ to 10 GHz
clear all;

y1 = logspace(8,10, 100);
noise_power = 10*log10((242 + 30)*y1*1.3806e-23);

figure(1);

plot (y1, noise_power, 'r')
title ('Noise Power vs bandwidth frequency');

xlabel('Bandwidth (Hz)') 
ylabel('Noise Power (dBW)') 

%% Noise Power vs bandwidth frequency 1MHz to 10 MHZ
clear all;
y1 = logspace(6,7,20);
noise_power = 10*log10((242 + 30)*y1*1.3806e-23);

figure(2);

plot (y1, noise_power, 'b')
title ('Noise Power vs bandwidth frequency 1MHz to 100 MHZ (242 K)');
grid on;

xlabel('Bandwidth (Hz)') 
ylabel('Noise Power (dBW)') 