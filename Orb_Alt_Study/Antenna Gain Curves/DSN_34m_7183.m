function [DS_HGA_X] = DSN_34m_7183()
% This function creates a spline curve of the dB gain versus angle off bore
% axis for the steered HGA on the Curiosity rover

G_bor = [0 68.41;
    0.1 68;
    0.2 68];

G_bor(:,1) = deg2rad(G_bor(:,1));
DS_HGA_X.lim = G_bor(end,1);
DS_HGA_X.bor = fit(G_bor(:,1),G_bor(:,2),'poly2');
