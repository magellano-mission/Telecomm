
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

%% NOTES
% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% SHARED INPUTS
%Antenna gain curves vs boresight angle and angle limit information
[UHF_G] = MRO_UHF_401_6();   %[dBi] MRO / Curiosity UHF antenna curve
%[LGA_X] = MRO_LGA_7183();    %[dBi] MRL X band LGA  curve
%[MGA_X] = MSL_MGA_7183();    %[dBi] MSL X band MGA curve
%[HGA_X] = Cur_HGA_7183();    %[dBi] Curiosity X band HGA curve
%[S_HGA_X] = MRO_HGA_7183();  %[dBi] MRO X band HGA curve
%[DSN_MGA_X] = DSN_34m_7183(); %[dBi] DSN 34m dish gain (curve not available)

%Environment parameters
T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space (REVISIT)
a_Earth = astroConstants(2);
a_Mars  = 1.52*a_Earth;

%Electra specifications for UHF baseline
elec.powr = [-170 -100];    %[dBW] Acceptable signal power range to feed Electra
elec.M = 4;                 %[bit/sym] Max modulation efficiency (BPS,QPS,etc)
elec.R = 0.5;               %[-] Error coding efficiency
elec.P = 'circ';            %[-]Polarisation of signal, 'circ' or 'lin'
elec.F = 3.9;                         %[dB] Noise Factor for Electra (half-duplex)
elec.Ts = T_amb*(10^(elec.F/10)-1);   %[K] System effective noise temperature

%SDST specifications for X and Ka band baseline
sdst.powr = [-188 -100];     %[dBW] Acceptable signal power range to feed SDST
sdst.M = 8;                  %[bit/sym] Max modulation efficiency (BPS,QPS,etc)
sdst.R = 0.5;                %[-] Error coding efficiency
sdst.P = 'circ';             %[-]Polarisation of signal, 'circ' or 'lin'
sdst.F = 2.5;                         %[dB] Noise Factor for SDST
sdst.Ts = T_amb*(10^(sdst.F/10)-1);   %[K] System effective noise temperature

%DSN specifications for X and Ka band baseline
dsn.powr = [-300 -50];
dsn.M = 8;
dsn.R = 0.5;
dsn.P = 'circ';
dsn.Ts = 20;


%% INPUTS - TYPICAL ORBITAL USER UHF STUDY
%Ground user is kept at fixed equatorial location with fixed hardware and 
%fixed power.Eclipses are not considered. The orbital altitude of the
%orbiter is varied
altitudes = 100:25:3000;              %[km] altitude range of interest
a = altitudes + astroConstants(24);     %[km] semi-major axis range of interest
keps = zeros(length(a),6);              %[km & rads] 
keps(:,1) = 1*a';                       %[km & rads]

ustat = [300,0,75,0,0,0];             %[km alt & degrees] typical user position
ustat(1) = ustat(1) + astroConstants(24);
ustat(2:6) = deg2rad(ustat(2:6));     %degrees to radians

dt = 10;                %[s] time step
sols = 3;               
t = 0: dt : sols*88620; %[s] Mars day = 88620, Earth day = 86400

sit = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter 
                  % to Mars orbiter, 3 - Mars to Earth (generic)
frq = 401.6e6;    %[Hz] carrier signal frequency
BW   = 1e6;       %[Hz] Bandwidth
powt = 40;        %[W] ground user RF power emitted
point = [0 0];    %[tran recv] 1 - antenna steered, 0 - not steered

%% FUNCTION - TYPICAL ORBITAL USER UHF STUDY
tic
out = struct;
res = zeros(length(a),7);
for i = 1:length(a)
[out(i).all,tots] = pass_over(1,sit,1,frq,BW,powt,elec,keps(i,:),ustat,t,dt,point,...
                  UHF_G,UHF_G,'Orbiter to Orbiter UHF',0);
    %storing totals 
    res(i,:) = [altitudes(i) tots(1) tots(2) tots(3) tots(4) tots(5) sols]; 
                            %[kJ]    [uJ]    [Gb]    vis     con      
end
toc

%% SAVE RESULTS
%save results and variables to output file for post-processing convenience
save('Output Files\UHF_comms_trends_sit2_out','out','res','ustat','keps',...
     'dt','sols','sit','frq','BW','powt','point','t')

%% POST-PROCESSING
%Load information from file if required 
%load('Output Files\UHF_comms_trends_sit2_out')
plots = zeros(length(out),6);
plots(:,1) = res(:,4)./sols;             %[Gb/sol] Daily data transfer
plots(:,2) = res(:,4)./res(:,2)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,3) = res(:,5)./sols;             %[hrs/sol] Visible Time
plots(:,4) = res(:,6)./sols;             %[hrs/sol] Downlink Time


%% STUDY PLOTTING
figure(1)
subplot(2,2,1)
plot(res(:,1),plots(:,1)*1000)
xlabel('Altitude [km]')
ylabel('Data Transfer [Mb/sol]')
title({'Orbiter Altitude vs Daily Data Transfer'})

subplot(2,2,2)
plot(res(:,1),plots(:,2)*1000)
xlabel('Altitude [km]')
ylabel('Data Transfer Efficiency [Mb/kJ/sol]')
title({'Orbiter Altitude vs Data Transfer Efficiency'})

subplot(2,2,3)
plot(res(:,1),plots(:,3))
xlabel('Altitude [km]')
ylabel('Visible Time  [hrs/sol]')
title({'Orbiter Altitude vs Daily Visible Time'})

subplot(2,2,4)
plot(res(:,1),plots(:,4))
xlabel('Altitude [km]')
ylabel('Downlink Time  [hrs/sol]')
title({'Orbiter Altitude vs Daily Downlink Time'})

sgtitle({'Unsteered UHF Comms w/ Electra','i=75 Orbital User to i=0 Equatorial Orbiter','40W & 2dBi Boresight Gain'})
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

return

%% PARTICULAR CASE PLOTTING 
%choose an altitude of interest to produce plots for
alt_int = 200;     %[km] altitude
[log,row] = ismember(alt_int,altitudes);

%choose a new time frame to produce the plots for
dt_plot = 10;
sols_plot = 1;
t_plot = 0:dt_plot:88620*sols_plot;

if log == 1
[a,~] = pass_over(1,sit,2,frq,BW,powt,elec,keps(row,:),ustat,t_plot,dt_plot,point,...
                  UHF_G,UHF_G,'Detailed Behaviour',1);
end