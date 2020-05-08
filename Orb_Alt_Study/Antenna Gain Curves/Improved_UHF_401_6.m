function [UHF] = Improved_UHF_401_6(plotting)
% This function creates spline curves of the dB gain versus angle off bore
% axis and elevation angle for the MRO and Curiosity UHF antennas 
% transmitting at 401.6 MHz

G_bor = [0 5;
    5 4.9;
    10 4.6;
    15 4.2;
    20 3.5;
    25 3.0;
    30 2.3;
    35 1.3;
    40 0.3;
    45 -0.7
    50 -2.0;
    55 -3.5;
    60 -4.9;
    65 -6.4;
    70 -8.0];

if plotting == 1
    UHF.plotting = fit(G_bor(:,1),G_bor(:,2),'smoothingspline');
    plot(UHF.plotting)
    xlabel('Boresight Offset (degrees)')
    ylabel('Directional Gain [dBi]')
end

G_bor(:,1) = deg2rad(G_bor(:,1));
UHF.lim = G_bor(end,1);
UHF.bor = fit(G_bor(:,1),G_bor(:,2),'smoothingspline');