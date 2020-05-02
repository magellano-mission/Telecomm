%% PROBLEM SET UP
clearvars
clc
clf

mars_day = 88620;

%% SHARED INPUTS
T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space
BW    = 1e6;       %[Hz] Assuming 1MHz for now
powt = 10;         %[W] ground user RF power emitted

%Electra specifications - using as a baseline
powr = [-170 -70]; %[dBW] Acceptable power range to feed Electra
code.M = [2 4];  %[bits/symbol] Signal level range (BPSK,QPSK) aka Modulation efficiency
code.R = [0.5 1];  %[-] Error coding efficiency
F = 3.9;           %[dB] Noise Factor for Electra (half-duplex)
T_s = T_amb*(10^(F/10)-1);  %System effective noise temperature

%Orbital altitudes and time span information
alts = [200, 600, 1000]; %[km] orbiter altitudes to analyse
dt = 30;                              %[s] time step
t = 0: dt : 0.2*mars_day;           %[s] time span being analysed

%% CASE 1 - MRO/Curiosity UHF at 401.6 MHz
frq_1 = 401.6e6;                       %[Hz] carrier signal frequency
[~, UHF_G_ele] = MRO_UHF_401_6();    %[dBi] antenna gain vs elevation

figure(1)
[UHF] = pass_over_power(frq_1,powt,powr,alts,t,dt,UHF_G_ele,'MRO/Curiosity LGA UHF',T_s,BW,code);

%% CASE 2 - MRO LGA at 7183 MHz
frq_2 = 7183e6;                       %[Hz] carrier signal frequency
[~, LGA_G_ele] = MRO_LGA_7183();    %[dBi] antenna gain vs elevation

figure(2)
[LGA] = pass_over_power(frq_2,powt,powr,alts,t,dt,LGA_G_ele,'MRO LGA X-Band',T_s,BW,code);