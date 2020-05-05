clearvars; clc
format compact

%% NOTES
% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% SHARED INPUTS
%Antenna gain curves vs boresight angle and angle limit information
[UHF_G] = MRO_UHF_401_6();   %[dBi] MRO / Curiosity UHF antenna curve
[LGA_X] = MRO_LGA_7183();    %[dBi] MRL X band LGA  curve
[MGA_X] = MSL_MGA_7183();    %[dBi] MSL X band MGA curve
[HGA_X] = Cur_HGA_7183();    %[dBi] Curiosity X band HGA curve

T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space (REVISIT)

dt = 10;                     %[s] time step
mars_day = 88620;            %[s] length of Mars day
t = 0: dt : 1*mars_day;      %[s] time span being analysed

%Electra specifications for UHF baseline
elec.powr = [-170 -100];    %[dBW] Acceptable signal power range to feed Electra
elec.M = 4;                 %[bit/sym] Max modulation efficiency (BPS,QPS,etc)
elec.R = 0.5;               %[-] Error coding efficiency
elec.P = 'circ';            %[-]Polarisation of signal, 'circ' or 'lin'
elec.F = 3.9;                         %[dB] Noise Factor for Electra (half-duplex)
elec.Ts = T_amb*(10^(elec.F/10)-1);   %[K] System effective noise temperature

%SDST specifications for X and Ka band baseline
sdst.powr = [-188 -100];     %[dBW] Acceptable signal power range to feed SDST
sdst.M = 16;                 %[bit/sym] Max modulation efficiency (BPS,QPS,etc)
sdst.R = 0.5;                %[-] Error coding efficiency
sdst.P = 'circ';             %[-]Polarisation of signal, 'circ' or 'lin'
sdst.F = 2.5;                         %[dB] Noise Factor for SDST
sdst.Ts = T_amb*(10^(sdst.F/10)-1);   %[K] System effective noise temperature    

%Orbiter PSEUDO-Keplerian elements (ALTITUDE, NOT SEMI-MAJOR AXIS)
%Create a constellation by entering the initial Keplerian elements of all orbiters
keps = [2000,0,15,0,0,0;      %[km,deg,deg,deg,deg,deg]
        2000,0,15,90,0,0;     %[a,e,i,OM,om,th,mu]
        2000,0,15,180,0,0;
        2000,0,15,270,0,0;
        17032,0,0,0,0,0];
keps(:,2:6) = deg2rad(keps(:,2:6));           %degrees to radians
keps(:,1)   = keps(:,1) + astroConstants(24); %altitude to radius in [km]

%Ground user initial position
lint = [10,0];               %[deg lat, deg long]
lint = deg2rad(lint);       %degrees to radians

%% CASE 1 - MRO/Curiosity UHF at 401.6 MHz
frq_1 = 401.6e6;    %[Hz] carrier signal frequency
powt_1 = 10;        %[W] ground user RF power emitted
BW_1   = 1e6;       %[Hz] Bandwidth
point = [0 0];      %[tran recv] 1 if antenna steered, 0 if not
[UHF] = pass_over(1,1,frq_1,powt_1,elec,keps,lint,t,dt,point,UHF_G,UHF_G,'MRO/Curiosity UHF',BW_1);

%% CASE 2 - Curiosity HGA sending to MSL MGA at 7183 MHz
frq_2 = 7183e6;     %[Hz] carrier signal frequency
powt_2 = 10;        %[W] ground user RF power emitted
BW_2 = 1e6;         %[Hz] Bandwidth
point = [1 0];      %[tran recv] 1 if antenna steered, 0 if not
[HGA] = pass_over(2,1,frq_2,powt_2,sdst,keps,lint,t,dt,point,HGA_X,MGA_X,'Curiosity HGA to MSL MGA X-Band',BW_2);

%% CASE 3 - Something to do with Phased array antennas
