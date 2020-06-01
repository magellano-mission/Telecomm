function [recv,tr,res] = pass_over(sit,frq,powt,hard,rstat,tstat,t,dt,plots)
%% PURPOSE
%This function takes inputs about the transmission method, hardware,
%spacecraft, land vehicles and planets being considered and reports plots
%and total values for some criteria.

%% INPUTS
%sit.cas    - [-]       1 - Mars ground to orbiter, 2 - Mars orbiter to
                    %orbiter, 3 - Mars to Earth
%frq    - [Hz]      Carrier wave frequency
%powt   - [W]       Transmitted RF power
%hard.powr   - [W,W]     Minimum and maximum signal power limits of receiver
%rstat   - [km,rad]  Keplerian elements of receiver 
%tstat  - [rad,rad] OR [km,rad] Initial latitude and longitude or keplerian elements of transmitter
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
if sit.cas == 1     %Transmitter on Mars surface
    %finding states over time vector
    [tr.st] = state_prop('ground',tstat,t,sit.cas,hard.tiltt);  %[km & km/s]
end

if sit.cas == 2     %Transmitter in Mars orbit
    %finding states over time vector
    [tr.st] = state_prop('orbiter',tstat,t,sit.cas,hard.tiltt);  %[km & km/s]
    tr.st.pnt = hard.dirt * tr.st.pnt;
end

if sit.cas == 3     %Transmitting from Mars to Earth
    %finding states of Mars over the time vector
    [tr.st] = state_prop('orbiter',tstat,t,sit.cas,hard.tiltt);  %[km & km/s]
end
    
%% ORBITAL USER STATES AND LINK CALCULATIONS
%Intialising structures and counters
recv = struct;
sumEner = 0; sumData = 0; 
cnt = -1;        %counter for intervals when a good connection is possible
vis = -1;        %counter for when the user is visible to the system
test1 = 0;
test2 = 0;

for i = 1:size(rstat,1)
    %Retrieving state vector and variable history
    [recv(i).st] = state_prop('orbiter',rstat(i,:),t,sit.cas,hard.tiltr); %[km & km/s]
     recv(i).st.rel = recv(i).st.pos - tr.st.pos;       %[m] receiver to transmitter relative position vector
     recv(i).st.rng = sqrt(sum(recv(i).st.rel'.^2))';   %[m] range between receiver and transmitter
     recv(i).st.rangrat = gradient(recv(i).st.rng,t);   %[m/s] range rate between receiver and transmitter
     recv(i).st.dopfrq = frq.*(1-recv(i).st.rangrat./physconst('LightSpeed'))-frq;  %[Hz] Doppler shift
     recv(i).st.pnt = hard.dirr * recv(i).st.pnt;
     
     %Initialising blank vectors
     recv(i).st.ele  = NaN(length(t),1);      %local elevation
     recv(i).st.dE   = zeros(length(t),1);    %time step energy transmission
     recv(i).st.PdB  = NaN(length(t),1);      %received signal power
     recv(i).st.CNR  = NaN(length(t),1);      %CNR
     recv(i).st.EbNo = NaN(length(t),1);      %Bit energy to noise ratio
     recv(i).st.datrat = zeros(length(t),1);  %time step data rate
     recv(i).st.data = zeros(length(t),1);    %time step data volume
     recv(i).st.pow = zeros(length(t),1);     %time step operating power

     for j=1:length(t)
            switch sit.cas
                case 1  %surface transmitter to orbital receiver
                    % find receiver elevation from relative distance triangle
                    [C] = elevate(recv(i).st.rng(j),tr.st.dis(j),recv(i).st.dis(j));
                    % check if receiver sat is above surface user horizon
                    if C > pi/2
                        recv(i).st.ele(j) = C - pi/2; %[rad] elevation angle
                        vis = vis + 1;                %transmitter is visible to receiver 
                        test1 = 1;
                    else % receiver is not visible
                        recv(i).st.rng(j) = NaN;      %remove value
                        recv(i).st.rangrat(j) = NaN;  %remove value
                        recv(i).st.dopfrq(j) = NaN;   %remove value
                    end
                    
                case  2  %orbital transmitter to orbital receiver
                    % find receiver elevation from relative distance triangle
                    [C] = elevate(recv(i).st.rng(j),tr.st.dis(j),recv(i).st.dis(j));
                    recv(i).st.ele(j) = C - pi/2;
                    %ele_block = -asin(astroConstants(24)*1000/tr.st.dis(j)); % apex half angle
                    % find distance when blocking by Mars occurs
                    dis_block = sqrt(recv(i).st.dis(j)^2 - (astroConstants(24)*1000)^2) ...
                        + sqrt(tr.st.dis(j)^2 - (astroConstants(24)*1000)^2);
                    
                    % check if range is compatible with visibility
                    if recv(i).st.rng(j) < dis_block
                    %if recv(i).st.ele(j) > ele_block
                        vis = vis + 1;
                        test1 = 1;
                    else  % satellites can't see each other
                       recv(i).st.ele(j) = NaN;
                       recv(i).st.rng(j) = NaN;
                       recv(i).st.rangrat(j) = NaN;
                       recv(i).st.dopfrq(j) = NaN;
                    end
                    
                case 3  % assuming Earth always visible to ECS
                    recv(i).st.ele(j) = pi/2;
                    test1 = 1;
                    test2 = 1;
            end
            
            % If visibility conditions are satisfied
            if test1 == 1
                gt_el = NaN; gr_el = NaN;
                
                % if antennas are pointed towards each other   
                if test2 == 0 && dot(tr.st.pnt(j,:),recv(i).st.pnt(j,:)) <= 0 || sit.cas == 3
            
                    %Find gain of transmission antenna
                    switch hard.point(1)
                        case 0          %transmission antenna is unpointed
                            borang = atan2(norm(cross(recv(i).st.rel(j,:),tr.st.pnt(j,:))),...
                                                  dot(recv(i).st.rel(j,:),tr.st.pnt(j,:)));
                            if borang >= hard.limst(1) && borang <= hard.limst(2)
                                gt_el = hard.gaint.bor(borang);    %[dBi]
                            else
                                gt_el = NaN;
                            end
                            
                        case 1          %transmission antenna is pointed
                            gt_el = hard.gaint.bor(0);
                    end
            
                    %Find gain of receiver antenna
                    switch hard.point(2)
                        case 0          %receiver antenna is unpointed
                                        %relative position vector is reversed
                            borang = atan2(norm(cross(-recv(i).st.rel(j,:),recv(i).st.pnt(j,:))),...
                                                  dot(-recv(i).st.rel(j,:),recv(i).st.pnt(j,:)));
                            if borang >= hard.limsr(1) && borang <= hard.limsr(2)
                                gr_el = hard.gainr.bor(borang);    %[dBi]
                            else
                                gr_el = NaN;
                            end
                        case 1          %receiver antenna is pointed
                            gr_el = hard.gainr.bor(0);
                        case 3          %ECS receiver is an unpointed prism
                            % calculate the angle between ECS and RS in
                            % abosolute, non-rotating coordinates
                            ang = wrapTo2Pi(atan2(recv(i).st.rel(j,2),recv(i).st.rel(j,1)));
                            % retrieve the gain from the 2D assembly gain
                            % profile
                            gr_el = interp1(hard.prism(:,1),hard.prism(:,2),ang);
                    end
                end
                if isnan(gt_el) == 0 && isnan(gr_el) == 0
                    %Calculate the energy and data transmitted between elements
                        [out] = zfun_link_rate(frq,recv(i).st.rng(j),recv(i).st.ele(j),powt,hard.cont,gt_el,gr_el,dt,sit);
                        recv(i).st.dE(j) = out.dE;        %[J] received energy   
                        recv(i).st.PdB(j) = out.PdB;      %[dBW] received power
                        recv(i).st.CNR(j) = out.CNRdB;
                        recv(i).st.EbNo(j) = out.ENdB;
                        recv(i).st.datrat(j) = out.datrat;
                        recv(i).st.data(j) = out.data;
                        recv(i).st.pow(j) = powt*dt;
                        if recv(i).st.datrat(j) > 0
                            cnt = cnt + 1;
                        end
                end
               
            end
     end
     recv(i).st.cumE = cumsum(recv(i).st.dE);         %[J] cumulative RF energy received per orbiter
     sumEner = sumEner + recv(i).st.cumE(end);     %[J] RF energy received by the whole constellation
     recv(i).st.cumDat = cumsum(recv(i).st.data);     %[bits] cumulative data received per orbiter
     sumData = sumData + recv(i).st.cumDat(end);   %[bits] cumulative data received by whole constellation
     recv(i).st.cumpow = cumsum(recv(i).st.pow);      %[J] cumulative RF power transmitted
end
     

%% totals - CUMULATIVE OUTPUTS
res(1) = (cnt.*dt.*powt)/1e3;       %[kJ] total energy consumed
res(2) = sumEner*1e6;               %[uJ] total RF energy received
res(3) = sumData/1e9;               %[Gb] total data transferred
res(4) = vis.*dt / 3600;            %[hrs] total visible time
res(5) = cnt.*dt / 3600;            %[hrs] total downlink time

%% PLOTS ('plots' = 1 if plots are required)
if plots == 1

switch sit.cas
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
     plot(time,rad2deg(recv(i).st.ele))
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Elevation [deg]')
     ylim([-90 90])
     grid on
end
hold off

subplot(3,3,2)
hold on
for i = 1:length(recv)
     plot(time,recv(i).st.rng/1000)
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
     plot(time,recv(i).st.rangrat/1000)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Range Rate [km/s]')
     grid on
end
hold off

subplot(3,3,4)
hold on
for i = 1:length(recv)
     plot(time,recv(i).st.dopfrq/1000)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Doppler Shift [kHz]')
     grid on
end
hold off

subplot(3,3,5)
hold on
for i = 1:length(recv)
     plot(time,recv(i).st.cumpow/1000)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Cumulative Transmitted Energy [kJ]')
     grid on
end
    
subplot(3,3,6)
hold on
for i = 1:length(recv)
     plot(time,recv(i).st.PdB)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Signal Power [dBW]')
     ylim(hard.cont.powr)
     grid on
end
hold off

if sit.cas == 3
   subplot(3,3,7)
   hold on
   for i = 1:length(recv)
     plot(time,recv(i).st.EbNo)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Eb/No [dB]')
     ylim([-2 inf])
     grid on
   end
else
    subplot(3,3,7)
    hold on
    for i = 1:length(recv)
        plot(time,recv(i).st.CNR)
        xlabel(tlabel)
        xlim([0 time(end)])
        ylabel('CNR [dB]')
        ylim([0 inf])
        grid on
    end
end
hold off

subplot(3,3,8)
hold on
for i = 1:length(recv)
     plot(time,recv(i).st.datrat/1e6)
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
     plot(time,recv(i).st.cumDat/1e9)
     xlabel(tlabel)
     xlim([0 time(end)])
     ylabel('Cumulative Data [Gb]')
     grid on
end
hold off
end
