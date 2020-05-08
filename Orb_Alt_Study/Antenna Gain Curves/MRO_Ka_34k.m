function [DS_HGA_Ka] = MRO_Ka_34k(plotting)
% This function creates a spline curve of the dB gain versus angle off bore
% axis for the steered HGA on the Curiosity rover

G_bor = [0 56.4;
    0.1 52.4;
    0.2 36.4;
    0.3 42.4;
    0.4 32.5;
    0.5 32.0];

if plotting == 1
    HGA.plot = fit(G_bor(:,1),G_bor(:,2),'smoothingspline');
    plot(HGA.plot)
    xlabel('Boresight Offset (degrees)')
    ylabel('Directional Gain [dBi]')
end


G_bor(:,1) = deg2rad(G_bor(:,1));
DS_HGA_Ka.lim = G_bor(end,1);
DS_HGA_Ka.bor = fit(G_bor(:,1),G_bor(:,2),'smoothingspline');
