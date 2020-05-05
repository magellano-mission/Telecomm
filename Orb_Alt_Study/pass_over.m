function [recv] = pass_over(cas,go,frq,powt,hard,keps,lint,t,dt,point,gaint,gainr,title,BW)

% Change nomenclature to range and range-rate to reduce confusion when
% navigation is introduced

%% INPUTS
%cas    - [-]       Case number
%go     - [-]       1 - run the case, 0 - don't run the case
%frq    - [Hz]      Carrier wave frequency
%powt   - [W]       Transmitted RF power
%hard.powr   - [W,W]     Minimum and maximum signal power limits of receiver
%keps   - [km,rad]  Keplerian elements of orbits to be studied 
%lint   - [rad,rad] Initial latitude and longitude of ground user (MAY ADD VELOCITY IN THE FUTURE)
%t      - [s]       Vector of time span to be analysed
%dt     - [s]       Time step
%gaint  - [dBi,rad] Gain vs angle curve of transmission antenna
%gainr  - [dBi,rad] Gain vs angle curve of transmission antenna
%title  - [str]     Title for plot window
%hard.Ts- [K]       Effective noise temperature of the receiving system
%BW     - [Hz]      Receiver bandwidth (MAY MAKE THIS DYNAMIC IN FUTURE)
%hard.M - [-]       Max modulation efficiency
%hard.R - [-]       Error encoding rate (0.5 for turbo coding)
%hard.P - [str]     Polarisation of receiver and signal 'lin' or 'cir'

%% GO / NO GO
if go == 0
    recv = NaN;
    return   
end

%% TRANSMITTER STATES

%finding state for ground user over the time vector
[tran] = state_prop('ground',lint,t);  %[km & km/s]
tran.pos = tran.pos*1000;               %[m] convert units
tran.vel = tran.vel*1000;               %[m/s] convert units                    
%ground user absolute distance from Mars centre
tran.dis = sqrt(sum(tran.pos'.^2))';    %[m]

%% ORBITAL USER STATES AND LINK CALCULATIONS
recv = struct;     %structure for storing orbiter state info

%finding state for the different altitude orbiters over the time vector
for i = 1:size(keps,1)
    %Retrieving state vector history and the half-chord length with Mars surface
    [out] = state_prop('orbiter',keps(i,:),t); %[km & km/s]
    recv(i).pos = out.pos*1000;         %[m,m,m] convert units
    recv(i).vel = out.vel*1000;         %[m/s,m/s,m/s] convert units
    recv(i).hchord = out.hchord*1000;   %[m] convert units
    
    %Relative distance in non-rotating, Mars centre reference frame
     recv(i).rel = recv(i).pos - tran.pos;         %[m] orbiter to ground  user relative position vector
     recv(i).dis = sqrt(sum(recv(i).pos'.^2))';   %[m] orbiter ro Mars Centre distance
     recv(i).rng = sqrt(sum(recv(i).rel'.^2))';   %[m] range from orbiter to ground user
     recv(i).rangrat = gradient(recv(i).rng);     %[m/s] range rate from orbiter to ground user
     
     %Initialising blank vectors
     recv(i).ele  = NaN(length(t),1);      %local elevation
     recv(i).dE   = zeros(length(t),1);    %time step energy transmission
     recv(i).PdB  = NaN(length(t),1);      %received signal power
     recv(i).CNR  = NaN(length(t),1);      %CNR
     recv(i).datrat = zeros(length(t),1);  %time step data transmission
     
     for j=1:length(t)
        if recv(i).rng(j) < recv(i).hchord    %if the range between user and satellite is compatible with visibility
            
            %identifying variables of the triangle then solving for angle C
            a = recv(i).rng(j);      %Range, orbiter to ground user
            b = tran.dis(j);         %Mars centre to ground user distance
            c = recv(i).dis(j);      %Mars centre to orbiter distance
            C = acos((a^2 + b^2 -c^2)/(2*a*b));
            %solving for elevation angle
            recv(i).ele(j) = C - pi/2;           %[rad]
            
            %finding single antenna gain for the given elevation angle
            if point(1) == 0    %transmission antenna is unpointed
                borang = abs(recv(i).ele(j) - pi/2);
                gt_el = gaint.bor(borang);    %[dBi]
            end
            if point(1) == 1    %transmission antenna is pointed
                gt_el = gaint.bor(0);   
            end
            if point(2) == 0    %receiver antenna is unpointed
                borang = abs(recv(i).ele(j)-pi/2);
                gr_el = gainr.bor(borang);    %[dBi]
            end
            if point(2) == 1    %receiver antenna is pointed
                gr_el = gainr.bor(0);
            end
            %calculating the total energy transmission per time step
            [out] = zfun_link_rate(frq,a,powt,hard,gt_el,gr_el,BW,dt);
            recv(i).dE(j) = out.dE;
            recv(i).PdB(j) = out.PdB;
            recv(i).CNR(j) = out.CNRdB;
            recv(i).datrat(j) = out.datrat;
        else
            recv(i).rng(j) = NaN;    %if satellite is not visible, erase relative distance
            recv(i).rangrat(j) = NaN;
        end
     end
     recv(i).cumE = cumsum(recv(i).dE,1); %vector for cumulative energy trans
     recv(i).data = cumsum(recv(i).datrat,1);
end

%% PLOTS
figure(cas)
subplot(3,3,1)
hold on
for i = 1:length(recv)
     plot(t,rad2deg(recv(i).ele))
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Elevation [deg]')
     ylim([0 90])
     grid on
end
hold off

subplot(3,3,2)
hold on
for i = 1:length(recv)
     plot(t,recv(i).rng/1000)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Range [km]')
     grid on
end
hold off

subplot(3,3,3)
hold on
for i = 1:length(recv)
     plot(t,recv(i).rangrat/1000)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Range Rate [km/s]')
     grid on
end
hold off

subplot(3,3,4)
hold on
for i = 1:length(recv)
     plot(t,recv(i).dE)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Received Power [W]')
     grid on
end
hold off

subplot(3,3,5)
hold on
for i = 1:length(recv)
     plot(t,recv(i).cumE)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Cumulative Energy [J]')
     grid on
end
    
subplot(3,3,6)
hold on
for i = 1:length(recv)
     plot(t,recv(i).PdB)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Signal Power [dBW]')
     ylim(hard.powr)
     grid on
end

subplot(3,3,7)
hold on
for i = 1:length(recv)
     plot(t,recv(i).CNR)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('CNR [dB]')
     grid on
end

subplot(3,3,8)
hold on
for i = 1:length(recv)
     plot(t,recv(i).datrat/1e6)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Data Rate [Mb/s]')
     grid on
end
hold off

subplot(3,3,9)
hold on
for i = 1:length(recv)
     plot(t,recv(i).data/1e6)
     xlabel('Time [s]')
     xlim([t(1) t(end)])
     ylabel('Cumulative Data [Mb]')
     grid on
end
sgtitle(title)
hold off
