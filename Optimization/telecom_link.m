% COMMUNICATION UPLINK
clc, clear, close all

% X-Band definition
f = 8.44e9;         %[Hz] 
c = 3e8;            %[m/s]
lambda = c / f;     %[m]

% Path loss
R_mars = 3389.5e3;                                  %[m]
path_loss = @(x) 20 * log10(4 * pi * x / lambda);   %[dB]
SMA = linspace(100e3, 20000e3, 1000) + R_mars;      %[m]
losses = path_loss(SMA - R_mars);                   %[dB]
design_margin = 3;                                  %[dB]

% Path loss at 10500km
% SMA_selected = 5000e3;                      %[m]
% loss = path_loss(SMA_selected - R_mars);    %[dB]

% Noise power
k = 1.3806e-23;     %[?]
T_eq = 250;         %[K]

% Modulation of communication
mods = 2;

% Data rate to achieve:
data_rate = 2e6;                    %[bps]
bandwidth = data_rate ./ mods;      %[Hz]

% Minimum CNR to achieve:
CNRs = [10.6 13.6];

% Noise Power:
P_noise = 10 * log10(k * T_eq .* bandwidth);        %[dBW]

% Minimum Recieved Power:
P_rec = CNRs + P_noise;     %[dBW]

% Reciever Gain
G_sat = 19;            %[dBi]

% USER EIRP Needed:
EIRP = zeros(2,length(losses));

for kk = 1 : length(CNRs)
EIRP(kk, :) = P_rec(kk) - G_sat + losses + design_margin;  %[dBW]
end

% USER Power needed wrt gains:
G_user = [5, 10, 15];

bits = [1 2];
figure
for jj = 1 : length(CNRs)
    subplot(1,2,jj)
    for kk = 1 : 3
        P_sat = 10.^((EIRP(jj,:) - G_user(kk))./10);
        plot(SMA/1000, P_sat, 'LineWidth', 2)
        hold on
    end
    % ylim([0, 50])
    xlim([3000, 12000])
    xline(R_mars/1e3);
    xline(10500);
    title(strcat(num2str(bits(jj)), ' Mbps'))
    xlabel('SMA [km]')
    ylabel('Power Emitted [W]')
    legend('G = 5', 'G = 10', 'G = 15')
    grid on
end

figure
plot(SMA/1e3, EIRP,  'LineWidth', 2)
xlim([3000, 12000])
xline(R_mars/1e3);
xline(10500);
title('User EIRP needed @8.44GHz with 1Mhz bandwidth')
xlabel('SMA [km]')
ylabel('EIRP [dBW]')
legend('1 Mbps (BPSK)', '2 Mbps (QPSK)')
grid on

%% COMMUNICATION DOWNLINK
clc, clear, close all

% X-Band definition
f = 8.5e9;         %[Hz] 
c = 3e8;            %[m/s]
lambda = c / f;     %[m]

% Path loss
R_mars = 3389.5e3;                                  %[m]
path_loss = @(x) 20 * log10(4 * pi * x / lambda);   %[dB]
SMA = linspace(100e3, 10000e3, 1000) + R_mars;      %[m]
losses = path_loss(SMA - R_mars);                   %[dB]
design_margin = 3;                                  %[dB]

% Path loss at 10500km
% SMA_selected = 5000e3;                      %[m]
% loss = path_loss(SMA_selected - R_mars);    %[dB]

% Noise power
k = 1.3806e-23;     %[?]
T_eq = 250;         %[K]

% Modulation of communication
mods = [2];

% Data rate to achieve:
data_rate = [8e6];                    %[bps]
bandwidth = data_rate ./ mods;      %[Hz]

% Minimum CNR to achieve:
CNRs = 10.6;

% Noise Power:
P_noise = 10 * log10(k * T_eq .* bandwidth);        %[dBW]

% Minimum Recieved Power:
P_rec = CNRs + P_noise;     %[dBW]

% Reciever Gain
G_user = 5;            %[dBi]


% USER EIRP Needed:
EIRP = zeros(2,length(losses));

for kk = 1 : length(data_rate)
EIRP(kk, :) = P_rec(kk) - G_user + losses + design_margin;  %[dBW]
end

% USER Power needed wrt gains:
G_sat = linspace(8, 40, 100);

bw = sqrt(30000 ./ 10.^(G_sat./10));

%%
% figure
% hold on
% for jj = 1 : length(data_rate)
%     for kk = 1 : length(G_sat)
%         P_sat = 10.^((EIRP(jj,:) - G_sat(kk))./10);
%         plot(SMA/1000, P_sat, 'LineWidth', 2)
%     end
% end

P_sat = zeros(length(losses), length(G_sat));

for jj = 1 : length(data_rate)
    for kk = 1 : length(G_sat)
        P_sat(:,kk) = 10.^((EIRP(jj,:) - G_sat(kk))./10);
    end
end

SMa = SMA./1000;
surf(SMa, bw, P_sat')

%%
% ylim([0, 50])
xlim([3000, 12000])
xline(R_mars/1e3);
xline(10500);
title('User power needed @8.44GHz with 1Mhz bandwidth')
xlabel('SMA [km]')
ylabel('Power Emitted [W]')
legend('bw = 97', 'bw = 54', 'bw = 17', 'bw = 5.4')
grid on

figure
plot(SMA/1e3, EIRP(1,:),  'LineWidth', 2)
xlim([3000, 12000])
xline(R_mars/1e3);
xline(10500);
title('User EIRP needed @8.44GHz with 1Mhz bandwidth')
xlabel('SMA [km]')
ylabel('EIRP [dBW]')
grid on


%% COMMUNICATION DOWNLINK
clc, clear, close all

% X-Band definition
f = 4.01e8;         %[Hz] 
c = 3e8;            %[m/s]
lambda = c / f;     %[m]

% Path loss
R_mars = 3389.5e3;                                  %[m]
path_loss = @(x) 20 * log10(4 * pi * x / lambda);   %[dB]
SMA = linspace(100e3, 10000e3, 1000) + R_mars;      %[m]
losses = path_loss(SMA - R_mars);                   %[dB]
design_margin = 3;                                  %[dB]

% Path loss at 10500km
% SMA_selected = 5000e3;                      %[m]
% loss = path_loss(SMA_selected - R_mars);    %[dB]

% Noise power
k = 1.3806e-23;     %[?]
T_eq = 250;         %[K]

% Modulation of communication
mods = [2];

% Data rate to achieve:
data_rate = [2e6];                    %[bps]
bandwidth = data_rate ./ mods;      %[Hz]

% Minimum CNR to achieve:
CNRs = 10.6;

% Noise Power:
P_noise = 10 * log10(k * T_eq .* bandwidth);        %[dBW]

% Minimum Recieved Power:
P_rec = CNRs + P_noise;     %[dBW]

% Reciever Gain
G_user = 5;            %[dBi]


% USER EIRP Needed:
EIRP = zeros(2,length(losses));

for kk = 1 : length(data_rate)
EIRP(kk, :) = P_rec(kk) - G_user + losses + design_margin;  %[dBW]
end

% USER Power needed wrt gains:
G_sat = linspace(8, 40, 100);

bw = sqrt(30000 ./ 10.^(G_sat./10));

%
figure
hold on
for jj = 1 : length(data_rate)
    for kk = 1 : length(G_sat)
        P_sat = 10.^((EIRP(jj,:) - G_sat(kk))./10);
        plot(SMA/1000, P_sat, 'LineWidth', 2)
    end
end
% ylim([0, 50])
xlim([3000, 12000])
xline(R_mars/1e3);
xline(10500);
title('User power needed @8.44GHz with 1Mhz bandwidth')
xlabel('SMA [km]')
ylabel('Power Emitted [W]')
legend('bw = 97', 'bw = 54', 'bw = 17', 'bw = 5.4')
grid on

figure
plot(SMA/1e3, EIRP(1,:),  'LineWidth', 2)
xlim([3000, 12000])
xline(R_mars/1e3);
xline(10500);
title('User EIRP needed @8.44GHz with 1Mhz bandwidth')
xlabel('SMA [km]')
ylabel('EIRP [dBW]')
grid on