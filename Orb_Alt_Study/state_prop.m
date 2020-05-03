function [pos, vel, ang, hchord] = state_prop(type,state,time)
%Function to propagate the states of surface and orbital vehicles over a
% given time period. A circular orbit is assumed for orbital vehicles.
% Functionality may be added later to allow surface user velocity.

% States are populated using a NON-ROTATING MARS CENTRED REFERENCE FRAME
% Rotations are defined positive counter-clockwise

%% INPUTS
%type      - [str] 'ground' or 'orbiter'
%state     - [rad & km] initial latitude and longitude for ground users, or
             %vector of initial keplerian elements for orbiters
%time      - [s] time vector over which to evaluate

vel = NaN;
hchord = NaN;

r_mars = astroConstants(24);      %[km] radius at surface of Mars
mu     = astroConstants(14);      %[km^3/s^2] Mars gravitational parameter
mars_day = 88620;                 %[s]
  
%% Surface User
if ismember("ground",type) == 1
    n = 2*pi / mars_day;          %[rad/s] Mars angular velocity about axis
    lat = state(1);               %[rad] user latitude
    lon = state(2);               %[rad] user longitude
    
    %initialise position and velocity matrices
    spp = zeros(length(time),3);  %[rad,rad,km] spherical position vector  [azi,ele,radius] = [lat,lon,radius]  
    %spa = zeros(length(time),3);  %[rad/s2,rad/s2,km/s2] spherical acceleration vector
    pos = zeros(length(time),3);  %[km,km,km] Cartesian position vector
    vel = zeros(length(time),3);  %[km/s,km/s,km/s] Cartesian velocity vector
    
    %propagate position vector
    spp(1,:) = [lat, lon, r_mars];          %user initial position
    spp(:,1) = n.*time' + spp(1,1);         %propagate azimuth
    spp(:,2) = spp(1,2);                    %propagate elevation (constant)
    spp(:,3) = r_mars;                      %propagate radius (constant)
    %spherical coords to Cartesian, angles +ve CCW from +ve x axis (right)
   [pos(:,1),pos(:,2),pos(:,3)] = sph2cart(spp(:,1),spp(:,2),spp(:,3));
    
    %propagate velocity vector
    spv = [n;0;0];                     %[rad/s,rad/s,km/s] spherical velocity vector
    speed = n*r_mars*cos(lat);         %instantaneous speed
    
    for i = 1:length(time) %I couldn't quickly think of a better way to do this
        dot = speed*sph2cartvec(spv,rad2deg(spp(i,1) + pi/2),0);
        vel(i,1) = dot(1);
        vel(i,2) = dot(2);
        vel(i,3) = 0;
    end
    
    ang = spp(:,1);
    
end
    
%% Orbital User
if alt > 0
    r = r_mars + alt;             %[km] orbital radius
    a = r;                        %[km] assuming circular orbit
    T = 2*pi*sqrt(a^3/mu);        %[s] orbital period
    n = 2*pi/T;                   %[rad/s] mean motion parameter
    
    spp = zeros(length(time),3);  %[km,rad,rad] spherical position vector
    pos = zeros(length(time),3);  %[km,km,km] cartesian position vector
    
    %propagate position vector
    spp(:,1) = n.*time';       %[rad] azimuth (offset from surface user)
    spp(:,2) = 0;                 %[km] elevation angle
    spp(:,3) = r;                 %[km] radius (constant)
    [pos(:,1),pos(:,2),pos(:,3)] = sph2cart(spp(:,1),spp(:,2),spp(:,3));
    pos = pos * 1000;             %[m] convert units to metres
    
    theta = 2*acos(r_mars/r);   %[rad] chord central angle  
    chord = 2*r*sin(theta/2);    %[km] chord length
    hchord = chord/2*1000;       %[m] convert units to metres and halve
    
    ang = spp(:,1);
  
end

