function [pnt_cart] = pointer(pos_spher,tilt)

%% PURPOSE
% This function takes as inputs the spherical position vector of the
% user vehicle and the tilts applied to the antenna direction relative to
% the nadir and zenith directions, and outputs the pointing direction time
% history in Cartesian coordinates

%% INPUTS
% pos_spher - [rad,rad,km] position vector history in spherical coordinates
% tilt      - [rad,rad]    tilts to be applied in the lat/elevation and
%                          long/azimuth directions

pnt = ones(length(pos_spher));
pnt(:,1) = pos_spher(:,1) + tilt(1);
pnt(:,2) = pos_spher(:,2) + tilt(2);

[pnt_cart(1,:),pnt_cart(2,:),pnt_cart(3,:)] = sph2cart(pnt(:,1),pnt(:,2),pnt(:,3));