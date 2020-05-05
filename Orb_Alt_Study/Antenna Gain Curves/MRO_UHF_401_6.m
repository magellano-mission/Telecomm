function [UHF] = MRO_UHF_401_6()
% This function creates spline curves of the dB gain versus angle off bore
% axis and elevation angle for the MRO and Curiosity UHF antennas 
% transmitting at 401.6 MHz

G_bor = [0 2.1;
    5 2.2;
    10 2.2;
    15 2.2;
    20 2.2;
    25 2.3;
    30 2.4;
    35 2.1;
    40 1.6;
    45 1.2
    50 0.5;
    55 -0.5;
    60 -1.0;
    65 -2.0;
    70 -2.8];

G_bor(:,1) = deg2rad(G_bor(:,1));
UHF.lim = G_bor(end,1);
UHF.bor = fit(G_bor(:,1),G_bor(:,2),'smoothingspline');



    
    