
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')

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


%% INPUTS - TYPICAL GROUND USER UHF STUDY
%Ground user is kept at fixed equatorial location with fixed hardware and 
%fixed power.Eclipses are not considered. The orbital altitude of the
%orbiter is varied
altitudes = 200:25:3000;              %[km] altitude range of interest
a = altitudes + astroConstants(24);     %[km] semi-major axis range of interest
keps = zeros(length(a),6);              %[km & rads] 
keps(:,1) = 1*a';                       %[km & rads]

ustat = [0,0];              %[deg lat, deg long] typical user position
ustat = deg2rad(ustat);     %degrees to radians

dt = 30;                %[s] time step
sols = 5;               
t = 0: dt : sols*88620; %[s] Mars day = 88620, Earth day = 86400

sit = 1;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter 
                  % to Mars orbiter, 3 - Mars to Earth (generic)
frq = 401.6e6;    %[Hz] carrier signal frequency
BW   = 1e6;       %[Hz] Bandwidth
powt = 15;        %[W] ground user RF power emitted
point = [0 0];    %[tran recv] 1 - antenna steered, 0 - not steered

%% FUNCTION
tic
out = struct;
res = zeros(length(a),6);
for i = 1:length(a)
[out(i).all,tots] = pass_over(1,sit,1,frq,BW,powt,elec,keps(i,:),ustat,t,dt,point,...
                  UHF_G,UHF_G,'Curiosity UHF to MRO UHF',0);
    %storing totals 
    res(i,:) = [altitudes(i) tots(1) tots(2) tots(3) tots(4) tots(5)]; 
                            %[kJ]    [uJ]    [Gb]    vis     con      
end
toc

%% POST-PROCESSING
plots = zeros(length(a),6);
plots(:,1) = res(:,4)./sols;             %[Gb/sol] Daily data transfer
plots(:,2) = res(:,4)./res(:,2)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,3) = res(:,5)./sols;             %[hrs/sol] Visible Time
plots(:,4) = res(:,6)./sols;             %[hrs/sol] Downlink Time


%% PLOTTING
subplot(2,2,1)
plot(res(:,1),plots(:,1))
xlabel('Altitude [km]')
ylabel('Data Transfer [Gb/sol]')
title({'Orbiter Altitude vs Daily Data Transfer'})

subplot(2,2,2)
plot(res(:,1),plots(:,2))
xlabel('Altitude [km]')
ylabel('Data Transfer Efficiency [Gb/kJ/sol]')
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

sgtitle({'Unsteered UHF Comms w/ Electra','1 Equatorial Ground User to 1 Equatorial Orbiter','15W & 2dBi Boresight Gain'})
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

