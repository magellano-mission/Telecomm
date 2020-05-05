function P = getSatellitePower(x, P_user_max, f, mods_sat, mods_user, data_rate_sat, data_rate_user, N)
% P:                power needed for satellite transmitter [W]
%
%   x:  [1x3]       [SMA (km), G_user (dB), bw (deg)] 
%   P_user_max:     maximum power for user transmitter [W]
%   f:              frequency needed [Hz]
%   mods_sat:       satellite downlink modulation 1 for BPSK, 2 for QPSK) []
%   mods_user:      user uplink modulation (1 for BPSK, 2 for QPSK) []
%   data_rate_sat:  satellite downlink data-rate [bps]
%   data_rate_user: user uplink data-rate [bps]
%   N:              number of users [ ]

if nargin == 7
    N = 1;
end

% Semimajor axis:
SMA = x(1);            %[km]

% Reciever Gain:
G_user = x(2);            %[dBi]

% Satellite beamwidth:
bw = x(3);  %[deg]

% Wavelength definition
c = 3e8;            %[m/s]
lambda = c / f;     %[m]

% Path loss
R_mars = 3389.5e3;                                  %[m]
path_loss = @(x) 20 * log10(4 * pi * x / lambda);   %[dB]
losses = path_loss(SMA*1e3 - R_mars);               %[dB]
design_margin = 3;                                  %[dB]

% Noise power data
k = 1.3806e-23;     %[?]
T_eq = 250;         %[K]

% Data rate to achieve to the user:
bandwidth_user = data_rate_sat ./ mods_sat;      %[Hz]

% Minimum CNR to achieve:
if mods_sat == 1
    CNRs_sat = 10.6;
elseif mods_sat == 2
    CNRs_sat = 13.6;
end

% Noise Power:
P_noise = 10 * log10(k * T_eq .* bandwidth_user);        %[dBW]

% Minimum Recieved Power:
P_rec = CNRs_sat + P_noise;     %[dBW]

% Satellite EIRP Needed:
EIRP = P_rec - G_user + losses + design_margin;  %[dBW]

% Gain of the satellite computation wrt beamwidth:
G_sat = 10 * log10(30000 / bw^2);

% Power needed for the satellite:
P_sat = 10.^((EIRP - G_sat)./10);


% Data rate to achieve at the satellite:
bandwidth_sat = data_rate_user ./ mods_user * N;      %[Hz]

% Minimum CNR to achieve:
if mods_user == 1
    CNRs_user = 10.6;
elseif mods_user == 2
    CNRs_user = 13.6;
end

% Noise Power:
P_noise_sat = 10 * log10(k * T_eq .* bandwidth_sat);        %[dBW]

% Minimum Recieved Power:
P_rec_sat = CNRs_user + P_noise_sat;     %[dBW]

% User Power Needed:
P_user = P_rec_sat - G_sat - G_user + losses + design_margin;  %[dBW]

if P_user > P_user_max
    P = Inf;
else
    P = P_sat;
end
end
