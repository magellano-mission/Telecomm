%% Received power vs distance @ 8GHz, 20db antennas
clear all;
close all;
distance = linspace (1000000,12000000,20);
carrier_freq = linspace (500e6,8e9,10);

sender_gain = 10;
sender_power = 15;
sender_EIRP = sender_gain + sender_power ;
receiver_gain = 20;
%received_power = -20*log10((4*pi*distance*8e9)/3e8) + sender_EIRP + receiver_gain ;

%set bandwidth @5MHz
noise_power = 10*log10((242 + 30)*5e6*1.3806e-23);
SNR_vec = [];

for i= 1:length(carrier_freq)
    for j = 1:length(distance)
        received_power = -20*log10((4*pi*distance(j)*carrier_freq(i)/3e8)) + sender_EIRP + receiver_gain ;
        SNR = received_power - noise_power;
        SNR_vec (i,j) = [SNR];
    end 
    plot (distance, SNR_vec(i,:));
    hold on;
end  

 
grid on;
plot (distance, SNR, 'g')
title ('SNR  vs altitude @ 5MHz');

xlabel('Distance (km)') 
ylabel('SNR')

annotation('arrow', [.3 .2], [.8 .3] );
dim = [.4 .5 .8 .3];
str = 'decreasing freq. 500MHz - 8GHz';
annotation('textbox',dim,'String',str,'FitBoxToText','on');