function [orb] = pass_over_power(frq,powt,powr,keps,lint,t,dt,gaint,gainr,title,T_s,BW,code)

% Change nomenclature to range and range-rate to reduce confusion when
% navigation is introduced

%% INPUTS
%frq    - [Hz]      Carrier wave frequency
%powt   - [W]       Transmitted RF power
%powr   - [W,W]     Minimum and maximum signal power limits of receiver
%keps   - [km,rad]  Keplerian elements of orbits to be studied 
%lint   - [rad,rad] Initial latitude and longitude of ground user (MAY ADD VELOCITY IN THE FUTURE)
%t      - [s]       Vector of time span to be analysed
%dt     - [s]       Time step
%gaint  - [dBi,rad] Gain vs angle curve of transmission antenna
%gainr  - [dBi,rad] Gain vs angle curve of transmission antenna
%title  - [str]     Title for plot window
%T_s    - [K]       Effective noise temperature of the receiving system
%BW     - [Hz]      Receiver bandwidth (MAY MAKE THIS DYNAMIC IN FUTURE)
%code.M - [-]       Max modulation efficiency
%code.R - [-]       Error encoding rate (0.5 for turbo coding)
%code.P - [str]     Polarisation of receiver and signal 'lin' or 'cir'

%% GROUND USER STATES
% Later this can just be changed to generic user, first need to think more
% about things like 3D angles between orbiters, Mars LOS blocking etc


%finding state for ground user over the time vector
[grd.pos,~,grd.ang,~] = state_prop('ground',lint,t);  %[m,m/s]
%ground user absolute distance from Mars centre
grd.dis = sqrt(sum(grd.pos'.^2))';        %[m]

%CHANGE TO KM HERE AFTER STATE SCRIPT

%% ORBITAL USER STATES AND LINK CALCULATIONS
orb = struct;     %structure for storing orbiter state info

%finding state for the different altitude orbiters over the time vector
for i = 1:length(keps)
    %Retrieving state vector history and the half-chord length with Mars surface
    [orb(i).pos,~,orb(i).ang,orb(i).hchord] = ...
        state_prop('orbiter',keps(i),t); %[m,m/s]
    
    %RELATIVE position between orbiter and ground user
    % IN MARS INERTIAL COORDINATES
     orb(i).rel = orb(i).pos - grd.pos;             %[m,m,m]
     
     %calculating the magitude of the orbiter-ground user distance and
     %orbiter-Mars centre distance
     orb(i).dis = sqrt(sum(orb(i).pos'.^2))';       %[m]
     orb(i).dis_rel = sqrt(sum(orb(i).rel'.^2))';   %[m]
     
     %calculating the local elevation angle between the ground user and
     % satellite when it is visible to the user
     orb(i).ele = NaN(length(t),1); %vector for local elevation
     orb(i).dE   = zeros(length(t),1); %vector for time step energy trans
     orb(i).PdB = NaN(length(t),1); %vector for CNR
     orb(i).CNR = NaN(length(t),1); %vector for CNR
     orb(i).rate   = zeros(length(t),1); %vector for time step energy trans
     for j=1:length(t)
        if orb(i).dis_rel(j) < orb(i).hchord    %if the distance between user and satellite is compatible with visibility
            
            %identifying variables of the triangle then solving for angle C
            a = orb(i).dis_rel(j);  %Relative distance, orbiter to ground user
            b = grd.dis(j);         %Mars centre to ground user distance
            c = orb(i).dis(j);      %Mars centre to orbiter distance
            C = acos((a^2 + b^2 -c^2)/(2*a*b));
            %solving for elevation angle
            orb(i).ele(j) = C - pi/2;           %[rad]
            %finding single antenna gain for the given elevation angle
            gt_el = gaint(orb(i).ele(j));    %[dBi]
            gr_el = gainr(orb(i).ele(j));    %[dBi]
            %calculating the total energy transmission per time step
            [out] = zfun_link_rate(frq,a,powt,powr,gt_el,gr_el,T_s,BW,dt,code);
            orb(i).dE(j) = out.E;
            orb(i).PdB(j) = out.PdB;
            orb(i).CNR(j) = out.CNRdB;
            orb(i).rate(j) = out.rate;
        else
            orb(i).dis_rel(j) = NaN;    %if satellite is not visible, erase relative distance
        end
     end
     orb(i).cumE = cumsum(orb(i).dE,1); %vector for cumulative energy trans
     orb(i).data = cumsum(orb(i).rate,1);
end

%% PLOTS
subplot(3,3,1)
hold on
for i = 1:length(orb)
     plot(t,rad2deg(orb(i).ele))
     xlabel('Time [s]')
     ylabel('Elevation [deg]')
end
hold off

subplot(3,3,2)
hold on
for i = 1:length(orb)
     plot(t,orb(i).dis_rel)
     xlabel('Time [s]')
     ylabel('Link Distance [m]')
end
hold off

subplot(3,3,3)
hold on
for i = 1:length(orb)
     plot(t,orb(i).dE)
     xlabel('Time [s]')
     ylabel('Received Power [W]')
end
hold off

subplot(3,3,4)
hold on
for i = 1:length(orb)
     plot(t,orb(i).cumE)
     xlabel('Time [s]')
     ylabel('Cumulative Energy [J]')
end
     %legend('200 km', '400 km', '600 km', '800 km','Location','northwest')

subplot(3,3,5)
hold on
for i = 1:length(orb)
     plot(t,orb(i).PdB)
     xlabel('Time [s]')
     ylabel('Signal Power [dBW]')
     ylim(powr)
end

subplot(3,3,6)
hold on
for i = 1:length(orb)
     plot(t,orb(i).CNR)
     xlabel('Time [s]')
     ylabel('CNR [dB]')
end

subplot(3,3,7)
hold on
for i = 1:length(orb)
     plot(t,orb(i).rate)
     xlabel('Time [s]')
     ylabel('Data Rate [bit/s]')
end
hold off

subplot(3,3,8)
hold on
for i = 1:length(orb)
     plot(t,orb(i).data)
     xlabel('Time [s]')
     ylabel('Cumulative Data [bits]')
end
sgtitle(title)
hold off
