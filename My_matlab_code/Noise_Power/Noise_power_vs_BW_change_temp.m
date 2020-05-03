%% Noise power vs bandwidth at 10 MHZ to 70 MHz (242 K, 30db LNA)
clear all;

y1 = linspace (10e6,70e6,10e6);
noise_power = 10*log10((242 + 30)*y1*1.3806e-23);

figure(6);

plot (y1, noise_power, 'o')
title ('Noise Power vs bandwidth frequency 10 MHZ to 70 MHz 242K ');

xlabel('Bandwidth (Hz)') 
ylabel('Noise Power (dBW)')

%%


y1 = linspace (10e6,70e6,10e6);
noise_power = 10*log10((500 + 30)*y1*1.3806e-23);

figure(7);

plot (y1, noise_power, 'o')
title ('Noise Power vs bandwidth frequency 10 MHZ to 70 MHz 500K ');

xlabel('Bandwidth (Hz)') 
ylabel('Noise Power (dBW)')