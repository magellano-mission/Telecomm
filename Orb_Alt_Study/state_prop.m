function [pos, vel, ang, hchord] = state_prop(alt,time)
%Function to propagate the states of surface and orbital vehicles over a
% given time period. It is assumed that surface vehicles and orbital
% vehicles begin the time period separate by 180 degrees around the equator
% of Mars

% ASSUMES CIRCULAR ORBIT

% INPUTS KM, OUTPUTS M
vel = NaN;
hchord = NaN;

r_mars = astroConstants(24);      %[km] radius at surface of Mars
mu     = astroConstants(14);      %[km^3/s^2] Mars gravitational parameter
mars_day = 88620;                 %[s]
n = 2*pi / mars_day;              %[rad/s] Mars rotational angular velocity
  
%% Surface User
if alt == 0
    %initialise blank matrices
    sph = zeros(length(time),3);  %[km,rad,rad] spherical position vector
    pos = zeros(length(time),3);  %[km,km,km] cartesian position vector
    
    %propagate position vector
    sph(:,1) = n.*time';          %[rad] azimuth
    %pos_rad(:,2) = NA            %[km] elevation angle
    sph(:,3) = r_mars;            %[km] radius (constant)
    [pos(:,1),pos(:,2),pos(:,3)] = sph2cart(sph(:,1),sph(:,2),sph(:,3));
    pos = pos * 1000;             %[m] convert units to metres
    
    ang = sph(:,1);
    
end
    
%% Orbital User

if alt > 0
    r = r_mars + alt;             %[km] orbital radius
    a = r;                        %[km] assuming circular orbit
    T = 2*pi*sqrt(a^3/mu);        %[s] orbital period
    n = 2*pi/T;                   %[rad/s] mean motion parameter
    
    sph = zeros(length(time),3);  %[km,rad,rad] spherical position vector
    pos = zeros(length(time),3);  %[km,km,km] cartesian position vector
    
    %propagate position vector
    sph(:,1) = n.*time'-pi;       %[rad] azimuth (offset from surface user)
    %pos_rad(:,2) =               %[km] elevation angle
    sph(:,3) = r;                 %[km] radius (constant)
    [pos(:,1),pos(:,2),pos(:,3)] = sph2cart(sph(:,1),sph(:,2),sph(:,3));
    pos = pos * 1000;             %[m] convert units to metres
    
    theta = 2*acos(r_mars/r);   %[rad] chord central angle  
    chord = 2*r*sin(theta/2);    %[km] chord length
    hchord = chord/2*1000;       %[m] convert units to metres and halve
    
    ang = sph(:,1);
  
end

