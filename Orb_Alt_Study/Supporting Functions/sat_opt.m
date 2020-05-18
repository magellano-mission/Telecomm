function [bw, Nsat] = sat_opt(sma, lat)

% Input:        - semi-major axes [km] 
%               - latitude range centered on equator e.g. [-30 +30] [deg]
%
% Output:       - beamwidth [deg]
%               - number of satellites needed 

% Setup
RM = 3393;              % [km]
h_sat = sma - RM;
lat = lat(2)*pi/180;

% Consider as angular radius the inscribed square diagonal
theta = lat*sqrt(2); 

% Beamwidth computation
gamma = abs(atan(sin(theta)/(1 + h_sat/RM - cos(theta))));
bw = 2*gamma;
bw = bw*180/pi;

Nsat = ceil(pi/theta);

