function [pnt] = pointer(pos_spher,tilt)

%% PURPOSE
% This function takes as inputs the zenith or nadir vector and the vector
% of the tilts that need to be applied and outputs the history of the
% antenna pointing vector for a nominally unsteered antenna

%% INPUTS
% pos_spher - [rad,rad,km] position vector history in spherical coordinates
% tilt      - [rad,rad] tilts to be applied in the lat/elevation and
%             long/azimuth directions

