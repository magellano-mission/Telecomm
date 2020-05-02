clearvars
clc

altitude = linspace(200,4000,3);     %[km] orbiter altitudes to analyse
mars_day = 88620;                   %[s]
dt = 10;                            %[s]
t = 0: dt : 0.33*mars_day;          %[s] time span being analysed

orb = struct;     %structure for storing orbiter state info
grd = struct;     %structure for storing ground user state info

%finding state for ground user over the time vector
[grd.pos,~,grd.ang,~] = state_prop(0,t);  %[m,m/s]
%ground user absolute distance from Mars centre
grd.dis = sqrt(sum(grd.pos'.^2))';

%finding state for the different altitude orbiters over the time vector
for i = 1:length(altitude)
    %Retrieving state vector history and the mars surface half-chord length
    % for the particular altitude orbit
    [orb(i).pos,~,orb(i).ang,orb(i).hchord] = state_prop(altitude(i),t); %[m,m/s]
    
    %calculating the RELATIVE position between the orbiter and ground user
    % IN MARS INERTIAL COORDINATES
     orb(i).rel = orb(i).pos - grd.pos;
     
     %calculating the magitude of the orbiter-ground user distance
     orb(i).dis = sqrt(sum(orb(i).pos'.^2))';
     orb(i).dis_rel = sqrt(sum(orb(i).rel'.^2))';
     
     %calculating the local elevation angle between the ground user and
     % satellite when it is visible to the user
     orb(i).ele = NaN(length(t),1); %blank matrix for local relative positions
     for j=1:length(t)
        if orb(i).dis_rel(j) < orb(i).hchord    %if the distance between user and satellite is compatible with visibility
            %identifying variables of the triangle then solving for angle C
            a = orb(i).dis_rel(j);  %Relative distanve orbiter to ground user
            b = grd.dis(j);         %Mars centre to ground user distance
            c = orb(i).dis(j);      %Mars centre to orbiter distance
            C = acos((a^2 + b^2 -c^2)/(2*a*b));
            orb(i).ele(j) = C - pi/2;   %elevation angle
        else
            orb(i).dis_rel(j) = NaN;    %if satellite is not visible, erase relative distance
        end
     end
     
     
     
end

%% PLOTS
subplot(1,2,1)
hold on
for i = 1:length(orb)
     plot(t,orb(i).ele)
     xlabel('Time [s]')
     ylabel('Elevation [rad]')
end
hold off

subplot(1,2,2)
hold on
for i = 1:length(orb)
     plot(t,orb(i).dis_rel)
     xlabel('Time [s]')
     ylabel('Link Distance [m]')
end
hold off

