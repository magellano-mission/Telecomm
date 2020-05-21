function [curve] = basic_antenna(type,gain_peak,HPBW,plotting)
% This function creates a basic quadratic curve for the antenna gain vs
% angle-off-boresight.

%% INPUTS
%gain_peak - [dBi] Peak gain of the antenna w.r.t. an isotropic antenna
%HPBW      - [deg] Half Power Beam Width of the antenna

if isequal(type,'helical') || isequal(type,'horn') || isequal(type,'parabolic') || isequal(type,'planar array')
    G_bor = [-deg2rad(HPBW/2), gain_peak - 3;
            0, gain_peak; 
            deg2rad(HPBW/2), gain_peak - 3];

    curve.bor = fit(G_bor(:,1),G_bor(:,2),'poly2');
end

if isequal(type,'phased array')
    G_bor = linspace(0,pi/2,1000)';
    G_bor(:,2) = gain_peak + 10.*log10((cos(G_bor(:,1))).^2.4);
    for i=1:length(G_bor)
        if G_bor(i,1) > pi/3
            G_bor(i,2) = -1000;
        end
    end
    curve.bor = fit(G_bor(:,1),G_bor(:,2),'linearinterp');
end

if plotting == 1
    plot(curve.bor)
    xlabel('Boresight Offset (radians)')
    xlim([0 pi/2])
    ylabel('Directional Gain [dBi]')
    ylim([0 inf])
end


    