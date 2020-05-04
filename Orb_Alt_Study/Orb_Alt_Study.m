clearvars; clc; clf

%% SHARED INPUTS
T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space
BW    = 1e6;       %[Hz] Assuming 1MHz for now
powt = 10;         %[W] ground user RF power emitted
dt = 30;                            %[s] time step
mars_day = 88620;                   %[s] length of Mars day
t = 0: dt : 0.2*mars_day;           %[s] time span being analysed

%Electra specifications for baseline
powr = [-170 -70];          %[dBW] Acceptable power range to feed Electra
code.M = 4;                 %[bit/sym] Max modulation efficiency (BPS,QPS,etc)
code.R = 0.5;               %[-] Error coding efficiency
code.P = 'circ';            %[-]Polarisation of signal, 'circ' or 'lin'
F = 3.9;                    %[dB] Noise Factor for Electra (half-duplex)
T_s = T_amb*(10^(F/10)-1);  %[K] System effective noise temperature

%Orbiter PSEUDO-Keplerian elements (ALTITUDE, NOT SEMI-MAJOR AXIS)
keps = [200,0,0,0,0,0;      %[km,deg,deg,deg,deg,deg]
        200,0,15,0,0,0;
        200,0,30,0,0,0];
keps(:,2:6) = deg2rad(keps(:,2:6));           %degrees to radians
keps(:,1)   = keps(:,1) + astroConstants(24); %altitude to radius in [km]

%Ground user initial position
lint = [15,0];               %[deg lat, deg long]
lint = deg2rad(lint);       %degress to radians

%% CASE 1 - MRO/Curiosity UHF at 401.6 MHz
frq_1 = 401.6e6;                       %[Hz] carrier signal frequency
[~, UHF_G_ele] = MRO_UHF_401_6();    %[dBi] antenna gain vs elevation

figure(1)
[UHF] = pass_over_power(frq_1,powt,powr,keps,lint,t,dt,UHF_G_ele,UHF_G_ele,'MRO/Curiosity LGA UHF',T_s,BW,code);

%% CASE 2 - MRO LGA at 7183 MHz
frq_2 = 7183e6;                       %[Hz] carrier signal frequency
[~, LGA_G_ele] = MRO_LGA_7183();    %[dBi] antenna gain vs elevation

figure(2)
[LGA] = pass_over_power(frq_2,powt,powr,keps,lint,t,dt,LGA_G_ele,LGA_G_ele,'MRO LGA X-Band',T_s,BW,code);