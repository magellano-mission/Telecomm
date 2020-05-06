function [DS_HGA_X] = MRO_HGA_7183()
% This function creates a spline curve of the dB gain versus angle off bore
% axis for the steered HGA on the Curiosity rover

G_bor = [0 46.7;
    0.1 46.4;
    0.2 45.7;
    0.3 44.4;
    0.4 42.5;
    0.5 39.7];

G_bor(:,1) = deg2rad(G_bor(:,1));
DS_HGA_X.lim = G_bor(end,1);
DS_HGA_X.bor = fit(G_bor(:,1),G_bor(:,2),'poly2');
