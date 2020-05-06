function [recv,totals] = pass_over(sit,num,frq,powt,hard,keps,ustat,t,dt,point,gaint,gainr,title,BW,plots)
%% PURPOSE
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
    %ground user absolute distance from Mars centre
    tran.dis = sqrt(sum(tran.pos'.^2))';    %[m]
end

if sit == 3     %Transmitter 'is' Mars
    %finding state for Mars over the time vector
    [tran] = state_prop('orbiter',ustat,t,sit);  %[km & km/s]
    tran.pos = tran.pos*1000;               %[m] convert units
    tran.vel = tran.vel*1000;               %[m/s] convert units                    
    %ground user absolute distance from Sun centre
    tran.dis = sqrt(sum(tran.pos'.^2))';    %[m]
end
    

%% ORBITAL USER STATES AND LINK CALCULATIONS

%Initialising structs
recv = struct;  
totals = struct;

%Intialising counters
sumEner = 0; 
sumData = 0; 
cnt = 0;        %counter for intervals when a good connection is possible
dcnt = 0;       %counter for intervals when a good connection is not possible
vis = 0;        %counter for when the user is visible to the system

for i = 1:size(keps,1)
    %Retrieving state vector history and the half-chord length with Mars surface
    [out] = state_prop('orbiter',keps(i,:),t,sit); %[km & km/s]
    recv(i).pos = out.pos*1000;         %[m,m,m] convert units
    recv(i).vel = out.vel*1000;         %[m/s,m/s,m/s] convert units
    if sit == 1
        recv(i).hchord = out.hchord*1000;   %[m] convert units
    end
    
    %Relative distance in non-rotating, Mars centre reference frame
     recv(i).rel = recv(i).pos - tran.pos;         %[m] orbiter to ground  user relative position vector
     recv(i).dis = sqrt(sum(recv(i).pos'.^2))';   %[m] orbiter ro Mars Centre distance
     recv(i).rng = sqrt(sum(recv(i).rel'.^2))';   %[m] range from orbiter to ground user
     recv(i).rangrat = gradient(recv(i).rng,t);     %[m/s] range rate from orbiter to ground user
     recv(i).dopfrq = frq.*(1 - recv(i).rangrat./physconst('LightSpeed')) - frq;
     
     %Initialising blank vectors
     recv(i).ele  = NaN(length(t),1);      %local elevation
     recv(i).dE   = zeros(length(t),1);    %time step energy transmission
     recv(i).PdB  = NaN(length(t),1);      %received signal power
     recv(i).CNR  = NaN(length(t),1);      %CNR
     recv(i).datrat = zeros(length(t),1);  %time step data transmission

     for j=1:length(t)
         switch sit     %Visibility test
             case 1     %Check if range is less than half chord length
                 y = recv(i).rng(j) < recv(i).hchord;
             case 2     %Check if Mars is not blocking the link, add later
                 y = 1; 
             case 3     %Check if Sun is not blocking the link, add later
                 y = 1;
         end
         
         
        if y == 1    %if visibility requirements are satisfied
            vis = vis + 1;
            switch sit
                case 1                   
                    %identifying variables of the triangle then solving for angle C
                    a = recv(i).rng(j);      %Range, orbiter to ground user
                    b = tran.dis(j);         %Mars centre to ground user distance
                    c = recv(i).dis(j);      %Mars centre to orbiter distance
                    C = acos((a^2 + b^2 -c^2)/(2*a*b));
                    %solving for elevation angle
                    recv(i).ele(j) = C - pi/2;           %[rad]
                case 2
                    recv(i).ele(j) = pi/2;
                    a = recv(i).rng(j);
                case 3
                    recv(i).ele(j) = pi/2;
                    a = recv(i).rng(j);
            end
            
            
            %Find gain of transmission antenna
            switch point(1)
                case 0          %transmission antenna is unpointed
                    borang = abs(recv(i).ele(j) - pi/2);
                    gt_el = gaint.bor(borang);    %[dBi]
                case 1          %transmission antenna is pointed
                    gt_el = gaint.bor(0);
            end
            
            %Find gain of receiver antenna
            switch point(2)
                case 0          %receiver antenna is unpointed
                    borang = abs(recv(i).ele(j)-pi/2);
                    gr_el = gainr.bor(borang);    %[dBi]
                case 1          %receiver antenna is pointed
                    gr_el = gainr.bor(0);
            end
            
            %Calculate the energy and data transmitted between elements
            [out] = zfun_link_rate(frq,a,powt,hard,gt_el,gr_el,BW,dt,sit);
            recv(i).dE(j) = out.dE;
            recv(i).PdB(j) = out.PdB;
            recv(i).CNR(j) = out.CNRdB;
            recv(i).datrat(j) = out.datrat;
            if recv(i).datrat(j) > 0
                cnt = cnt + 1;
            else
                dcnt = dcnt + 1;
            end
        else    %Visibility is not available
            recv(i).rng(j) = NaN;
            recv(i).rangrat(j) = NaN;
            recv(i).dopfrq(j) = NaN;
        end
     end
     recv(i).cumE = cumsum(recv(i).dE,1); %cumulative RF energy transmitted per orbiter
     sumEner = sumEner + recv(i).cumE(end); %RF energy transmitted to the whole constellation
     recv(i).data = cumsum(recv(i).datrat,1); %cumulative data transmitted per orbiter
     sumData = sumData + recv(i).data(end); %cumulative data transmitted to whole constellation
end
Power = cnt.*dt.*powt;
cTime = cnt.*dt;
dcTime = dcnt.*dt;
visTime = vis.*dt;
Energy = sumEner;
Data = sumData;
totals.energy_consumption_kJ = Power/1e3;
totals.energy_transmitted_uJ = Energy*1e6;
totals.data_transmitted_Gb = Data/1e6;
totals.efficiency_GB_per_kJ = (Data/1e6)/(Power/1e3);
totals.study_time_hrs = t(end) / 3600;
totals.visible_time_hrs = visTime / 3600;
totals.downlink_connected_time_hrs = cTime / 3600;
totals.downlink_disconnected_time_hrs = dcTime / 3600;

%% PLOTS ('plots' = 1 if plots are required)
if plots == 1

figure(num)
subplot(3,3,1)
hold on
for i = 1:length(recv)
     plot(t/3600,rad2deg(recv(i).ele))
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('Elevation [deg]')
     ylim([0 90])
     grid on
end
hold off

subplot(3,3,2)
hold on
for i = 1:length(recv)
     plot(t/3600,recv(i).rng/1000)
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('Range [km]')
     ylim([0 inf])
     grid on
end
hold off

subplot(3,3,3)
hold on
for i = 1:length(recv)
     plot(t/3600,recv(i).rangrat/1000)
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('Range Rate [km/s]')
     grid on
end
hold off

subplot(3,3,4)
hold on
for i = 1:length(recv)
     plot(t/3600,recv(i).dopfrq/1000)
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('Doppler Shift [kHz]')
     grid on
end
hold off

%subplot(3,3,4)
%hold on
%for i = 1:length(recv)
%     plot(t,recv(i).dE)
%     xlabel('Time [s]')
%     xlim([t(1) t(end)])
%     ylabel('Received Power [W]')
%     grid on
%end
%hold off

subplot(3,3,5)
hold on
for i = 1:length(recv)
     plot(t/3600,recv(i).cumE)
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('Cumulative Energy [J]')
     grid on
end
    
subplot(3,3,6)
hold on
for i = 1:length(recv)
     plot(t/3600,recv(i).PdB)
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('Signal Power [dBW]')
     ylim(hard.powr)
     grid on
end

subplot(3,3,7)
hold on
for i = 1:length(recv)
     plot(t/3600,recv(i).CNR)
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('CNR [dB]')
     grid on
end

subplot(3,3,8)
hold on
for i = 1:length(recv)
     plot(t/3600,recv(i).datrat/1e6)
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('Data Rate [Mb/s]')
     grid on
end
hold off

subplot(3,3,9)
hold on
for i = 1:length(recv)
     plot(t/3600,recv(i).data/1e6)
     xlabel('Time [hrs]')
     %xlim([t(1) t(end)/3600])
     ylabel('Cumulative Data [Mb]')
     grid on
end
sgtitle(title)
hold off
end
