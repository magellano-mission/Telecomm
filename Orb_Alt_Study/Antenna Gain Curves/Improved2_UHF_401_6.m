function [UHF] = Improved2_UHF_401_6(plotting)
% This function creates spline curves of the dB gain versus angle off bore
% axis and elevation angle for the MRO and Curiosity UHF antennas 
% transmitting at 401.6 MHz

G_bor = [0 5;
    5 5;%4.9;
    10 5;%4.6;
    15 5;%4.2;
    20 5;%3.5;
    25 5;%3.0;
    30 5;%2.3;
    35 0;%1.3;
    40 0;%0.3;
    45 0;%-0.7
    50 0;%-2.0;
    55 0;%-3.5;
    60 0;%-4.9;
    65 0;%-6.4;
    70 0];%-8.0

G_bor(:,2) = G_bor(:,2) + 3;

if plotting == 1
    UHF.plotting = fit(G_bor(:,1),G_bor(:,2),'smoothingspline');
    plot(UHF.plotting)
    xlabel('Boresight Offset (degrees)')
    ylabel('Directional Gain [dBi]')
end

G_bor(:,1) = deg2rad(G_bor(:,1));
UHF.lim = G_bor(end,1);
UHF.bor = fit(G_bor(:,1),G_bor(:,2),'smoothingspline');