%% Path loss vs distance 
distance = 1000000:500:12000000;

path_loss1 = 20*log10((4*pi*distance*400e6)/3e8);

figure (4);
plot (distance, path_loss1, 'g')
title ('Path loss vs distance');

xlabel('Distance (km)') 
ylabel('Path loss (dBW)')

hold on;



%1 GHz

path_loss3 = 20*log10((4*pi*distance*1e9)/3e8);
plot (distance, path_loss3, 'c')

xlabel('Distance (km)') 
ylabel('Path loss (dBW)')

% 4GHZ
path_loss2 = 20*log10((4*pi*distance*400e7)/3e8);
plot (distance, path_loss2, 'y')


xlabel('Distance (km)') 
ylabel('Path loss (dBW)')
legend ('400MHz', '4GHz')


%16 GHz

path_loss4 = 20*log10((4*pi*distance*16e9)/3e8);
plot (distance, path_loss4, 'r')

xlabel('Distance (km)') 
ylabel('Path loss (dBW)')


%32 GHz
path_loss5 = 20*log10((4*pi*distance*32e9)/3e8);
plot (distance, path_loss5, 'black')

xlabel('Distance (km)') 
ylabel('Path loss (dBW)')
legend ('400MHz', '1GHz', '4GHz', '16GHz', '32GHz')
hold off;