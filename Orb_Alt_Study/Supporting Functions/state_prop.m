function [out] = state_prop(type,state,time,sit)
%% PURPOSE
%This function propagates the states of surface and orbital vehicles over a
%given time period. Circular orbits are assumed. A NON-ROTATING Mars centred 
%reference frame is used with rotations defined positive counter-clockwise.
%Functionality can be added to allow elliptical orbits and non-zero surface
%velocities.

%% INPUTS
%type      - [str] 'ground' or 'orbiter'
%state     - [rad & km] initial latitude and longitude for ground users, or
             %vector of initial keplerian elements for orbiters
%time      - [s] time vector over which to evaluate
%sit       - [-] 1 - Mars ground user to Mars orbiter, 2 - Mars orbiter to
                    %Mars orbiter, 3 - Mars to Earth
                    
%% CASE VARIABLES                   
switch sit
    case 1 
        r_mars = astroConstants(24);      %[km] radius at surface of Mars
        mu     = astroConstants(14);      %[km^3/s^2] Mars gravitational parameter
        mars_day = 88620;                 %[s] length of Mars day
    case 2
        r_mars = astroConstants(24);      %[km] radius at surface of Mars
        mu     = astroConstants(14);      %[km^3/s^2] Mars gravitational parameter
    case 3
        mu = astroConstants(4);           %[km^3/s^2] Sun gravitational parameter
end

  
%% Mars Surface User
if ismember("ground",type) == 1
    n = 2*pi / mars_day;          %[rad/s] Mars angular velocity about axis
    lat = state(1);               %[rad] user latitude
    lon = state(2);               %[rad] user longitude
    
    %initialise position and velocity matrices
    spp = zeros(length(time),3);  %[rad,rad,km] spherical position vector  [azi,ele,radius] = [lat,lon,radius]  
    spv = zeros(length(time),3);  %[rad/s,rad/s,km/s] spherical velocity vector
    pos = zeros(length(time),3);  %[km,km,km] Cartesian position vector
    vel = zeros(length(time),3);  %[km/s,km/s,km/s] Cartesian velocity vector
    
    %propagate position vector
    spp(1,:) = [lon, lat, r_mars];          %user initial position [az,ele,r] w.r.t non-rotating Mars centre
    spp(:,1) = n.*time' + spp(1,1);         %propagate azimuth
    spp(:,2) = spp(1,2);                    %propagate elevation (constant)
    spp(:,3) = r_mars;                      %propagate radius (constant)
    %spherical coords to Cartesian, angles +ve CCW from +ve x axis (right)
   [pos(:,1),pos(:,2),pos(:,3)] = sph2cart(spp(:,1),spp(:,2),spp(:,3));
    
    %propagate velocity vector
    spv(:,1) = spp(:,1) + pi/2;             %perpendicular to radius
    spv(:,2) = 0;                           %change in elevation
    spv(:,3) = n*r_mars*cos(lat);           %scalar speed
    %spherical coords to Cartesian, angles +ve CCW from +ve x axis (right)
   [vel(:,1),vel(:,2),vel(:,3)] = sph2cart(spv(:,1),spv(:,2),spv(:,3));
   
   %Outputs for ground user
   out.pos = pos;                   %[km] Cartesian position vector
   out.vel = vel;                   %[km] Cartesian velocity vector
   %out.ang = spp(:,1);              %[rad] Azimuth w.r.t Mars centre
    
end
    
%% Orbital User
if ismember("orbiter",type) == 1
    r = state(1);                 %[km] orbital radius
    a = r;                        %[km] assuming circular orbit
    T = 2*pi*sqrt(a^3/mu);        %[s] orbital period
    n = 2*pi/T;                   %[rad/s] mean motion parameter
    
    %initialise position and velocity matrices
    pos = zeros(length(time),3);
    vel = zeros(length(time),3);
    
    %propagate true anomaly
    keps(1,1:6) = state;
    keps = repmat(keps,length(time),1);
    keps(:,5) = n.*time' + keps(1,5);
    
    %convert Keplerian elements to Cartesian
    for i = 1:length(time)
        [p, v] = par2car(keps(i,1),keps(i,2),keps(i,3),keps(i,4),keps(i,5),keps(i,6),mu);
        pos(i,:) = p';
        vel(i,:) = v';
    end
    
    %Outputs for orbiter
    out.pos = pos;               %[km] Cartesian position vector
    out.vel = vel;               %[km/s] Cartesian velocity vector
    
    if sit == 1
        %calculate the half-chord length for checking orbiter visibility
        theta = 2*acos(r_mars/r);    %[rad] chord central angle  
        chord = 2*r*sin(theta/2);    %[km] chord length
        hchord = chord/2;            %[m] halve to get one-way distance
        out.hchord = hchord;         %[km] Half chord length
    end
    
    if sit == 2
        %calculate the maximum visibility angle w.r.t. Mars
        ang = asin(r_mars/r);
        out.ang = ang;
    end
end

