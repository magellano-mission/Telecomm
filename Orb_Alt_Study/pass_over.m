function [orb] = pass_over(cas,go,frq,powt,powr,keps,lint,t,dt,gaint,gainr,title,T_s,BW,code)

% Change nomenclature to range and range-rate to reduce confusion when
% navigation is introduced

%% INPUTS
%cas    - [-]       Case number
%go     - [-]       1 - run the case, 0 - don't run the case
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

%% GO / NO GO
if go == 0
    orb = NaN;
    return   
end

%% GROUND USER STATES
% Later this can just be changed to generic user, first need to think more
% about things like 3D angles between orbiters, Mars LOS blocking etc


%finding state for ground user over the time vector
[grd] = state_prop('ground',lint,t);  %[km & km/s]
grd.pos = grd.pos*1000;               %[m] convert units
grd.vel = grd.vel*1000;               %[m/s] convert units                    
%ground user absolute distance from Mars centre
grd.dis = sqrt(sum(grd.pos'.^2))';    %[m]

%CHANGE TO KM HERE AFTER STATE SCRIPT

%% ORBITAL USER STATES AND LINK CALCULATIONS
orb = struct;     %structure for storing orbiter state info

%finding state for the different altitude orbiters over the time vector
for i = 1:size(keps,1)
    %Retrieving state vector history and the half-chord length with Mars surface
    [out] = state_prop('orbiter',keps(i,:),t); %[km & km/s]
    orb(i).pos = out.pos*1000;         %[m] convert units
    orb(i).vel = out.vel*1000;         %[m/s] convert units
    orb(i).hchord = out.hchord*1000;   %[m] convert units
    
    %Relative distance in non-rotating, Mars centre reference frame
     orb(i).rel = orb(i).pos - grd.pos;         %[m] orbiter to ground  user relative position vector
     orb(i).dis = sqrt(sum(orb(i).pos'.^2))';   %[m] orbiter ro Mars Centre distance
     orb(i).rng = sqrt(sum(orb(i).rel'.^2))';   %[m] range from orbiter to ground user
     
     %Initialising blank vectors
     orb(i).ele  = NaN(length(t),1);      %local elevation
     orb(i).dE   = zeros(length(t),1);    %time step energy transmission
     orb(i).PdB  = NaN(length(t),1);      %received signal power
     orb(i).CNR  = NaN(length(t),1);      %CNR
     orb(i).datrat = zeros(length(t),1);  %time step data transmission
     orb(i).rangrat = zeros(length(t),1); %time step data transmission
     
     for j=1:length(t)
        if orb(i).rng(j) < orb(i).hchord    %if the range between user and satellite is compatible with visibility
            %calculate the range rate
            orb_rangrat = dot(orb(i).rel(j),orb(i).vel(j));
            grd_rangrat = dot(orb(i).rel(j),grd.vel(j));
            orb(i).rangrat(j) = orb_rangrat + grd_rangrat;
            
            %identifying variables of the triangle then solving for angle C
            a = orb(i).rng(j);      %Range, orbiter to ground user
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
            orb(i).dE(j) = out.dE;
            orb(i).PdB(j) = out.PdB;
            orb(i).CNR(j) = out.CNRdB;
            orb(i).datrat(j) = out.datrat;
        else
            orb(i).rng(j) = NaN;    %if satellite is not visible, erase relative distance
            orb(i).rangrat(j) = NaN;
        end
     end
     orb(i).cumE = cumsum(orb(i).dE,1); %vector for cumulative energy trans
     orb(i).data = cumsum(orb(i).datrat,1);
end

%% PLOTS
figure(cas)
subplot(3,3,1)
hold on
for i = 1:length(orb)
     plot(t,rad2deg(orb(i).ele))
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Elevation [deg]')
     ylim([0 90])
     grid on
end
hold off

subplot(3,3,2)
hold on
for i = 1:length(orb)
     plot(t,orb(i).rng/1000)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Range [km]')
     grid on
end
hold off

subplot(3,3,3)
hold on
for i = 1:length(orb)
     plot(t,orb(i).rangrat/1000)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Range Rate [km/s]')
     grid on
end
hold off

subplot(3,3,4)
hold on
for i = 1:length(orb)
     plot(t,orb(i).dE)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Received Power [W]')
     grid on
end
hold off

subplot(3,3,5)
hold on
for i = 1:length(orb)
     plot(t,orb(i).cumE)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Cumulative Energy [J]')
     grid on
end
    
subplot(3,3,6)
hold on
for i = 1:length(orb)
     plot(t,orb(i).PdB)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Signal Power [dBW]')
     ylim(powr)
     grid on
end

subplot(3,3,7)
hold on
for i = 1:length(orb)
     plot(t,orb(i).CNR)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('CNR [dB]')
     grid on
end

subplot(3,3,8)
hold on
for i = 1:length(orb)
     plot(t,orb(i).datrat/1e6)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Data Rate [Mb/s]')
     grid on
end
hold off

subplot(3,3,9)
hold on
for i = 1:length(orb)
     plot(t,orb(i).data/1e6)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Cumulative Data [Mb]')
     grid on
end
sgtitle(title)
hold off
