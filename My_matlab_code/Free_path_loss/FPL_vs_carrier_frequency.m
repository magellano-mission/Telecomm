%% Signal Power vs carrier frequency
y2 = logspace(6,8);
free_path_loss = -20*log10((4*pi*11000000*y2)/3e8);

figure(3);
plot (y2, free_path_loss, 'b')
title ('Free path loss vs carrier frequency at 11,000km');

xlabel('Frequency (Hz)') 
ylabel('Path Loss (dBW)') 