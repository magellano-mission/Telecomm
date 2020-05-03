%% Received power vs distance @ 401MHz
distance = 1000000:500:12000000;

sender_gain = 5;
sender_power = 13;
sender_EIRP = sender_gain + sender_power;
receiver_gain = 20;
received_power = -20*log10((4*pi*distance*401e6)/3e8) + sender_EIRP + receiver_gain ;

figure (5);
plot (distance, path_loss1, 'g')
title ('Received power vs distance @ 401MHz');

xlabel('Distance (km)') 
ylabel('Received power (dBW)')
