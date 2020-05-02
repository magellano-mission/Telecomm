function [bor, ele] = MRO_LGA_7183()
% This function creates a spline curve of the dB gain versus angle off bore
% axis for the MRO LGA operating at 7183 MHz

G_bor = [0 0;
    5 -0.15;
    10 -0.3;
    15 -0.8;
    20 -1.2;
    25 -1.9;
    30 -2.5;
    35 -3.5;
    40 -4.3;
    45 -5.7;
    50 -6.5;
    55 -7.6;
    60 -8.5;
    65 -10.0;
    70 -11.5;
    75 -12.7;
    80 -14.0;
    85 -15.7;
    90 -16.6;];

G_bor(:,2) = G_bor(:,2) + 17;
G_bor(:,1) = deg2rad(G_bor(:,1));
G_ele = G_bor;
G_ele(:,2) = flip(G_ele(:,2));
    
bor = fit(G_bor(:,1),G_bor(:,2),'poly2');
ele = fit(G_ele(:,1),G_ele(:,2),'poly2');


    
    