function [out] = state_prop(type,state,time,sit,tilt)
%% PURPOSE
% This function propagates the states of surface and orbital vehicles over 
% a given time period. A NON-ROTATING Mars centred reference frame is used 
% with rotations defined positive counter-clockwise. Functionality can be 
% added to allow non-zero surface velocities.

%% INPUTS
%type      - [str] 'ground' or 'orbiter'
%state     - [rad & km] initial latitude and longitude for ground users, or
             % vector of initial keplerian elements for orbiters
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
    %antenna nominal pointing direction history in Cartesian coords
    [out.pnt] = pointer(spp,tilt);

    %propagate velocity vector
    spv(:,1) = spp(:,1) + pi/2;             %perpendicular to radius
    spv(:,2) = 0;                           %change in elevation
    spv(:,3) = n*r_mars*cos(lat);           %scalar speed
    %spherical coords to Cartesian, angles +ve CCW from +ve x axis (right)
   [vel(:,1),vel(:,2),vel(:,3)] = sph2cart(spv(:,1),spv(:,2),spv(:,3));
   
   %Outputs for ground user
   out.pos = pos*1000;                   %[m] Cartesian position vector
   out.vel = vel*1000;                   %[m] Cartesian velocity vector
   out.dis = sqrt(sum(out.pos'.^2))';    %[m] distance from Mars centre
    
end
    
%% Orbital User
if ismember("orbiter",type) == 1
    r = state(1);                 %[km] orbital radius
    a = r;                        %[km] assuming circular orbit
    T = 2*pi*sqrt(a^3/mu);        %[s] orbital period
    n = 2*pi/T;                   %[rad/s] mean motion parameter
    
    %initialise position and velocity matrices
    spp = zeros(length(time),3);  %[rad,rad,km]
    pos = zeros(length(time),3);  %[km,km,km]
    vel = zeros(length(time),3);  %[km/s,km/s,km/s]
    
    %propagate true anomaly
    keps(1,1:6) = state;
    keps = repmat(keps,length(time),1);
    %if the orbit is circular
    if state(2) == 0   
        keps(:,6) = n.*time' + keps(1,6);   %[rad] true anomaly
        keps(:,6) = wrapTo2Pi(keps(:,6));
    %if the orbit is elliptic
    else              
        keps(:,7) = n.*time' + keps(1,6);   %[rad] mean anomaly
        keps(:,7) = wrapTo2Pi(keps(:,7));   %[rad] wrapping to one circle
        %Initial eccentric anomaly
        E0 = 2*atan(sqrt((1-state(2))/(1+state(2)))*tan(keps(1,7)/2));
        
        for i = 2:length(time)
            %Creating function relating mean and eccentric anomalies
            ecc = @(E) E - state(2)*sin(E) - (E0 - state(2)*sin(E0)) - (keps(i,7) - keps(1,7));
            %Solving for eccentric anomaly at the given time step
            keps(i,8) = fzero(ecc,[0 2*pi+0.1]);
            %Calculating the true anomaly at the given time step
            keps(i,6) = wrapTo2Pi(2*atan(tan(keps(i,8)/2)/(sqrt((1-state(2))/(1+state(2))))));
        end
    end
    %convert Keplerian elements to Cartesian
    for i = 1:length(time)
        [p, v] = par2car(keps(i,1),keps(i,2),keps(i,3),keps(i,4),keps(i,5),keps(i,6),mu);
        pos(i,:) = p';
        vel(i,:) = v';
    end
    
    % antenna nominal pointing direction history in Cartesian coords
    [spp(:,1),spp(:,2),spp(:,3)] = cart2sph(pos(:,1),pos(:,2),pos(:,3));
    [out.pnt] = pointer(spp,tilt);
    
    %Outputs for orbiter
    out.pos = pos*1000;               %[m] Cartesian position vector
    out.vel = vel*1000;               %[m/s] Cartesian velocity vector
    out.anm = keps(:,6);              %[rad] true anomaly of orbiter
    out.dis = sqrt(sum(out.pos'.^2))';%[m] orbiter distance from focus
end

