function [recv,res] = pass_over(sit,frq,powt,hard,keps,ustat,t,dt,plots)
%% PURPOSE
%This function takes inputs about the transmission method, hardware,
%spacecraft, land vehicles and planets being considered and reports plots
%and total values for some criteria.

%% INPUTS
%sit    - [-]       1 - Mars ground to orbiter, 2 - Mars orbiter to
                    %orbiter, 3 - Mars to Earth
%frq    - [Hz]      Carrier wave frequency
%powt   - [W]       Transmitted RF power
%hard.powr   - [W,W]     Minimum and maximum signal power limits of receiver
%keps   - [km,rad]  Keplerian elements of orbits to be studied 
%ustat  - [rad,rad] Initial latitude and longitude of ground user (MAY ADD VELOCITY IN THE FUTURE)
%       - [km,rad] OR Keplerian elements for orbital vehicles or planets
%t      - [s]       Vector of time span to be analysed
%dt     - [s]       Time step
%hard.gaint  - [dBi,rad] Gain vs angle curve of transmission antenna
%hard.gainr  - [dBi,rad] Gain vs angle curve of transmission antenna
%hard.Ts- [K]       Effective noise temperature of the receiving system
%hard.BW- [Hz]      Receiver bandwidth (MAY MAKE THIS DYNAMIC IN FUTURE)
%hard.M - [-]       Max modulation efficiency
%hard.R - [-]       Error encoding rate (0.5 for turbo coding)
%hard.P - [str]     Polarisation of receiver and signal 'lin' or 'cir'

%% TRANSMITTER STATES
if sit == 1     %Transmitter on Mars surface
    %finding state for Mars ground user over the time vector
    [tran] = state_prop('ground',ustat,t,sit);  %[km & km/s]
    tran.pos = tran.pos*1000;               %[m] convert units
    tran.vel = tran.vel*1000;               %[m/s] convert units                    
    %ground user absolute distance from Mars centre
    tran.dis = sqrt(sum(tran.pos'.^2))';    %[m]
end

if sit == 2     %Transmitter in Mars orbit
    %finding state for Mars orbital user over the time vector
    [tran] = state_prop('orbiter',ustat,t,sit);  %[km & km/s]
    tran.pos = tran.pos*1000;               %[m] convert units
    tran.vel = tran.vel*1000;               %[m/s] convert units                    
    %orbital user absolute distance from Mars centre
    tran.dis = sqrt(sum(tran.pos'.^2))';    %[m]
end

if sit == 3     %Transmitter 'is' Mars
    %finding state for Mars over the time vector
    [tran] = state_prop('orbiter',ustat,t,sit);  %[km & km/s]
    tran.pos = tran.pos*1000;               %[m] convert units
    tran.vel = tran.vel*1000;               %[m/s] convert units                    
    %planet absolute distance from Sun centre
    tran.dis = sqrt(sum(tran.pos'.^2))';    %[m]
end
    
%% ORBITAL USER STATES AND LINK CALCULATIONS

%Intialising structures and counters
recv = struct;
sumEner = 0; sumData = 0; 
cnt = -1;        %counter for intervals when a good connection is possible
vis = -1;        %counter for when the user is visible to the system

for i = 1:size(keps,1)
    %Retrieving state vector history and the half-chord length with Mars surface
    [out] = state_prop('orbiter',keps(i,:),t,sit); %[km & km/s]
    recv(i).pos = out.pos*1000;         %[m,m,m] convert units
    recv(i).vel = out.vel*1000;         %[m/s,m/s,m/s] convert units
    recv(i).dis = sqrt(sum(recv(i).pos'.^2))';   %[m] orbiter to Mars Centre distance
    recv(i).anm = out.anm;
    
    %Relative distance in non-rotating, Mars centre reference frame
     recv(i).rel = recv(i).pos - tran.pos;         %[m] orbiter to ground  user relative position vector
     recv(i).rng = sqrt(sum(recv(i).rel'.^2))';   %[m] range from orbiter to ground user
     recv(i).rangrat = gradient(recv(i).rng,t);     %[m/s] range rate from orbiter to ground user
     recv(i).dopfrq = frq.*(1 - recv(i).rangrat./physconst('LightSpeed')) - frq;
     
     %Initialising blank vectors
     recv(i).ele  = NaN(length(t),1);      %local elevation
     recv(i).dE   = zeros(length(t),1);    %time step energy transmission
     recv(i).PdB  = NaN(length(t),1);      %received signal power
     recv(i).CNR  = NaN(length(t),1);      %CNR
     recv(i).datrat = zeros(length(t),1);  %time step data transmission
     recv(i).data = zeros(length(t),1);

     for j=1:length(t)
            switch sit
                case 1   %identifying variables of triangle and solving for elevation angle                
                    a = recv(i).rng(j);      %Range, orbiter to ground user
                    b = tran.dis(j);         %Distance, Mars centre to ground user distance
                    c = recv(i).dis(j);      %Distance, Mars centre to orbiter distance
                    C = acos((a^2 + b^2 -c^2)/(2*a*b));
                    if C > pi/2
                        recv(i).ele(j) = C - pi/2; %[rad] convert elevation to an angle w.r.t bore axis
                        vis = vis + 1;
                    else
                        recv(i).rng(j) = NaN;
                        recv(i).rangrat(j) = NaN;
                        recv(i).dopfrq(j) = NaN;
                    end
                case  2   %identifying variables of triangle and solving for elevation angle                
                    a = recv(i).rng(j);      %Range, orbiter to ground user
                    b = tran.dis(j);         %Distance, Mars centre to ground user distance
                    c = recv(i).dis(j);      %Distance, Mars centre to orbiter distance
                    C = acos((a^2 + b^2 -c^2)/(2*a*b));
                    if C > pi/2
                        recv(i).ele(j) = C - pi/2; %[rad] convert elevation to an angle w.r.t bore axis
                        vis = vis + 1;
                    else
                        recv(i).rng(j) = NaN;
                        recv(i).rangrat(j) = NaN;
                        recv(i).dopfrq(j) = NaN;
                    end    
                case 3
                    recv(i).ele(j) = pi/2;
                    a = recv(i).rng(j);
            end
            
            % Mainly for case 2, checking that the higher orbiter is within
            % view of the upper/rear of the lower orbiter
            if recv(i).ele(j) > 0 && recv(i).ele(j) < pi
            
                %Find gain of transmission antenna
                switch hard.point(1)
                    case 0          %transmission antenna is unpointed
                        borang = abs(recv(i).ele(j) - pi/2);
                        gt_el = hard.gaint.bor(borang);    %[dBi]
                    case 1          %transmission antenna is pointed
                        gt_el = hard.gaint.bor(0);
                end
            
                %Find gain of receiver antenna
                switch hard.point(2)
                    case 0          %receiver antenna is unpointed
                        borang = abs(recv(i).ele(j)-pi/2);
                        gr_el = hard.gainr.bor(borang);    %[dBi]
                    case 1          %receiver antenna is pointed
                        gr_el = hard.gainr.bor(0);
                    case 3
                        gr_el = interp1(hard.prism(:,1),hard.prism(:,2),recv(i).anm(j));
                end
            
                %Calculate the energy and data transmitted between elements
                [out] = zfun_link_rate(frq,a,powt,hard.cont,gt_el,gr_el,dt,sit);
                recv(i).dE(j) = out.dE;
                recv(i).PdB(j) = out.PdB;
                recv(i).CNR(j) = out.CNRdB;
                recv(i).datrat(j) = out.datrat;
                recv(i).data(j) = out.data;
                if recv(i).datrat(j) > 0
                    cnt = cnt + 1;
                end 
            end
     end
     recv(i).cumE = cumsum(recv(i).dE);       %cumulative RF energy transmitted per orbiter
     sumEner = sumEner + recv(i).cumE(end);     %RF energy transmitted to the whole constellation
     recv(i).cumDat = cumsum(recv(i).data);   %cumulative data transmitted per orbiter
     sumData = sumData + recv(i).cumDat(end);   %cumulative data transmitted to whole constellation
end
     

%% totals - CUMULATIVE OUTPUTS
res(1) = (cnt.*dt.*powt)/1e3;       %[kJ] total energy consumed
res(2) = sumEner*1e6;               %[uJ] total RF energy received
res(3) = sumData/1e9;               %[Gb] total data transferred
res(4) = vis.*dt / 3600;            %[hrs] total visible time
res(5) = cnt.*dt / 3600;            %[hrs] total downlink time
%disp(title)
%disp(recv.totals)

%% PLOTS ('plots' = 1 if plots are required)
if plots == 1

switch sit
    case 1
        time = t/3600;
        tlabel = 'Time [hrs]';
    case 2
        time = t/3600;
        tlabel = 'Time [hrs]';
    case 3
        time = t/86400;
        tlabel = 'Time [days]';
end
    
subplot(3,3,1)
hold on
for i = 1:length(recv)
     plot(time,rad2deg(recv(i).ele))
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Elevation [deg]')
     ylim([0 90])
     grid on
end
hold off

subplot(3,3,2)
hold on
for i = 1:length(recv)
     plot(time,recv(i).rng/1000)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Range [km]')
     ylim([0 inf])
     grid on
end
hold off

subplot(3,3,3)
hold on
for i = 1:length(recv)
     plot(time,recv(i).rangrat/1000)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Range Rate [km/s]')
     grid on
end
hold off

subplot(3,3,4)
hold on
for i = 1:length(recv)
     plot(time,recv(i).dopfrq/1000)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Doppler Shift [kHz]')
     grid on
end
hold off

%subplot(3,3,4)
%hold on
%for i = 1:length(recv)
%     plot(t,recv(i).dE)
%     xlabel('Time [s]')
%     ylabel('Received Power [W]')
%     grid on
%end
%hold off

subplot(3,3,5)
hold on
for i = 1:length(recv)
     plot(time,recv(i).cumE)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Cumulative Energy [J]')
     grid on
end
    
subplot(3,3,6)
hold on
for i = 1:length(recv)
     plot(time,recv(i).PdB)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Signal Power [dBW]')
     ylim(hard.cont.powr)
     grid on
end

subplot(3,3,7)
hold on
for i = 1:length(recv)
     plot(time,recv(i).CNR)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('CNR [dB]')
     ylim([-2 inf])
     grid on
end

subplot(3,3,8)
hold on
for i = 1:length(recv)
     plot(time,recv(i).datrat/1e6)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Data Rate [Mb/s]')
     ylim([0 inf])
     grid on
end
hold off

subplot(3,3,9)
hold on
for i = 1:length(recv)
     plot(time,recv(i).cumDat/1e9)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Cumulative Data [Gb]')
     grid on
end
%legend('RS Orbiter 1','RS Orbiter 2','RS Orbiter 3','RS Orbiter 4','Location','northwest')
hold off
end
