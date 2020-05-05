function [HGA_X] = Cur_HGA_7183()
% This function creates a spline curve of the dB gain versus angle off bore
% axis for the steered HGA on the Curiosity rover

G_bor = [0 25.5;
    2 24.1;
    5 20.4];

G_bor(:,1) = deg2rad(G_bor(:,1));
HGA_X.lim = G_bor(end,1);
HGA_X.bor = fit(G_bor(:,1),G_bor(:,2),'poly2');
