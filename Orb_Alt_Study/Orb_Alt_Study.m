clearvars; clc; clf
format compact

%% SHARED INPUTS
T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space
BW    = 1e6;       %[Hz] Assuming 1MHz for now
powt = 10;         %[W] ground user RF power emitted
dt = 10;                            %[s] time step
mars_day = 88620;                   %[s] length of Mars day
t = 0: dt : 1*mars_day;           %[s] time span being analysed

%Electra specifications for baseline
powr = [-170 -70];          %[dBW] Acceptable power range to feed Electra
code.M = 8;                 %[bit/sym] Max modulation efficiency (BPS,QPS,etc)
code.R = 0.5;               %[-] Error coding efficiency
code.P = 'circ';            %[-]Polarisation of signal, 'circ' or 'lin'
F = 3.9;                    %[dB] Noise Factor for Electra (half-duplex)
T_s = T_amb*(10^(F/10)-1);  %[K] System effective noise temperature

%Orbiter PSEUDO-Keplerian elements (ALTITUDE, NOT SEMI-MAJOR AXIS)
%You can create a constellation by entering the initial Keplerian elements
% of all the orbiters in the constellation
keps = [2000,0,15,0,0,0;      %[km,deg,deg,deg,deg,deg]
        2000,0,15,90,0,0;     %[a,e,i,OM,om,th,mu]
        2000,0,15,180,0,0;
        2000,0,15,270,0,0];
keps(:,2:6) = deg2rad(keps(:,2:6));           %degrees to radians
keps(:,1)   = keps(:,1) + astroConstants(24); %altitude to radius in [km]

%Ground user initial position
lint = [0,0];               %[deg lat, deg long]
lint = deg2rad(lint);       %degrees to radians

%% CASE 1 - MRO/Curiosity UHF at 401.6 MHz
frq_1 = 401.6e6;                       %[Hz] carrier signal frequency
[~, UHF_G_ele] = MRO_UHF_401_6();    %[dBi] antenna gain vs elevation
[UHF] = pass_over(1,1,frq_1,powt,powr,keps,lint,t,dt,UHF_G_ele,UHF_G_ele,'MRO/Curiosity LGA UHF',T_s,BW,code);

%% CASE 2 - MRO LGA at 7183 MHz
frq_2 = 7183e6;                       %[Hz] carrier signal frequency
[~, LGA_G_ele] = MRO_LGA_7183();    %[dBi] antenna gain vs elevation
[LGA] = pass_over(2,0,frq_2,powt,powr,keps,lint,t,dt,LGA_G_ele,LGA_G_ele,'MRO LGA X-Band',T_s,BW,code);

%% CASE 3 - MRO Improved 'MGA' at 7183 MHz
frq_3 = 7183e6;                       %[Hz] carrier signal frequency
[~, MGA_G_ele] = MRO_MGA_7183();    %[dBi] antenna gain vs elevation
[MGA] = pass_over(2,1,frq_3,powt,powr,keps,lint,t,dt,MGA_G_ele,MGA_G_ele,'MRO MGA X-Band',T_s,BW,code);
