function [DS_HGA_Ka] = DSN_34m_34k(plotting)
% This function creates a spline curve of the dB gain versus angle off bore
% axis for the steered HGA on the Curiosity rover

G_bor = [0 68.41;
    0.1 68;
    0.2 68];

if plotting == 1
end
G_bor(:,2) = G_bor(:,2) + 6.3;
G_bor(:,1) = deg2rad(G_bor(:,1));
DS_HGA_Ka.lim = G_bor(end,1);
DS_HGA_Ka.bor = fit(G_bor(:,1),G_bor(:,2),'poly2');
