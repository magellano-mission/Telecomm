clear all;
y1 = linspace (1e6,10e6, 20);
SNR = linspace (2, 12, 10);

data_rate_vec = [];
for i = 1:length (SNR)
   for j = 1:length(y1)
data_rate = y1(j)*log10(2*(1 + SNR(i))); 
data_rate_vec (i,j)  = [data_rate];
   end
plot (y1, data_rate_vec (i,:))


hold on;

%legend ([date_rate_vec(1,:), data_rate_vec (10,:)], 'SNR 2','SNR 12'));
end 


xlabel('Bandwidth (Hz)') 
ylabel('Data Rate (bps)')
annotation('arrow',  [.6 .5], [.3 .5]);
dim = [.4 .3 .5 .3];
str = 'increasing SNR 2 - 12';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
yticks(0:1e6:15e6);

