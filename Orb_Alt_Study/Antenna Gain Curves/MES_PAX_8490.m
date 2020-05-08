function [PAX] = MES_PAX_8490(plotting)
% This function creates a spline curve of the dB gain versus angle off bore
% axis for the steered HGA on the Curiosity rover

G_bor = [0 27.8;
    15 27.4;
    30 26.2;
    45 24.0;
    60 21.0];

if plotting == 1
    PAX.plot = fit(G_bor(:,1),G_bor(:,2),'poly2');
    plot(PAX.plot)
    xlabel('Boresight Offset (degrees)')
    ylabel('Directional Gain [dBi]')
end


G_bor(:,1) = deg2rad(G_bor(:,1));
PAX.lim = G_bor(end,1);
PAX.bor = fit(G_bor(:,1),G_bor(:,2),'poly2');
