%% Received power vs distance @ 8GHz, 20db antennas
distance = 1000000:500:12000000;

sender_gain = 10;
sender_power = 15;
sender_EIRP = sender_gain + sender_power ;
receiver_gain = 20;
received_power = -20*log10((4*pi*distance*8e9)/3e8) + sender_EIRP + receiver_gain ;

figure (7);
plot (distance, received_power, 'g')
title ('Received power vs distance @ 8GHz');

xlabel('Distance (km)') 
ylabel('Received power (dBW)')