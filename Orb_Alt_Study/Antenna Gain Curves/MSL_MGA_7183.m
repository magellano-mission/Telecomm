function [MGA_X] = MSL_MGA_7183(plotting)
% This function creates a spline curve of the dB gain off the cruise stage
% X-band MGA on the MSL/Curiosity mission

G_bor = [0 19;
    5 18;
    10 16;
    15 13;
    20 8;
    25 5;
    30 3;
    35 0;
    40 -10];

if plotting == 1
    MGA.plot = fit(G_bor(:,1),G_bor(:,2),'smoothingspline');
    plot(MGA.plot)
    xlabel('Boresight Offset (degrees)')
    ylabel('Directional Gain [dBi]')
end

G_bor(:,1) = deg2rad(G_bor(:,1));
MGA_X.lim = G_bor(end,1);
MGA_X.bor = fit(G_bor(:,1),G_bor(:,2),'poly2');
