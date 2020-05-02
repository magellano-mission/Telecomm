function [UHF_Gs] = MRO_UHF_401_6()
% This function creates a spline curve of the dB gain versus angle off bore
% axis for the MRO and Curiosity UHF antennas transmitting at 401.6 MHz
G_pat = [0 2.1;
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
    70 -2.8;
    75 -3.6;
    80 -4.4;
    85 -5.2;
    90 -6.0;];
    
UHF_Gs = fit(G_pat(:,1),G_pat(:,2),'smoothingspline');
    
    