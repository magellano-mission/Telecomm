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
[LGA_X] = MRO_LGA_7183();    %[dBi] MRL X band LGA  curve
[MGA_X] = MSL_MGA_7183();    %[dBi] MSL X band MGA curve
[HGA_X] = Cur_HGA_7183();    %[dBi] Curiosity X band HGA curve
[S_HGA_X] = MRO_HGA_7183();  %[dBi] MRO X band HGA curve
[DSN_MGA_X] = DSN_34m_7183(); %[dBi] DSN 34m dish gain (curve not available)

%Environment parameters
T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space (REVISIT)
dt = 10;                  %[s] time step
t = 0: dt : 88620;      %[s] Mars day = 88620, Earth day = 86400
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


%% (PROOF) CASE 1 - MRO/Curiosity UHF at 401.6 MHz
ustat1 = [4.6,0];              %[deg lat, deg long]
ustat1 = deg2rad(ustat1);     %degrees to radians
keps1 = [3788.15,0,92.89,0,0,0];     %[altitude km,deg,deg,deg,deg,deg]
                             %[a,e,i,OM,om,th,mu]

keps1(:,2:6) = deg2rad(keps1(:,2:6));           %degrees to radians
%keps1(:,1)   = keps1(:,1) + astroConstants(24); %altitude to radius in [km]
sit1 = 1;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter 
                   % to Mars orbiter, 3 - Mars to Earth (generic)
cas1 = 1;          %[-] case number
frq1 = 401.6e6;    %[Hz] carrier signal frequency
powt1 = 8.5;        %[W] ground user RF power emitted
BW1   = 1e6;       %[Hz] Bandwidth
point1 = [0 1];    %[tran recv] 1 - antenna steered, 0 - not steered
%Run the function and output graphs and totals
[Cur_UHF] = pass_over(0,sit1,cas1,frq1,powt1,elec,keps1,ustat1,t,dt,point1,...
                  UHF_G,UHF_G,'Curiosity UHF to MRO UHF',BW1,1);

%% CASE 2 - Orbiter to Orbiter at 401.6 MHz
ustat2 = [200,0,0,0,0,0];
keps2 = [500,0,15,0,0,0;     %[altitude km,deg,deg,deg,deg,deg]
        500,0,15,90,0,0;     %[a,e,i,OM,om,th,mu]
        500,0,15,180,0,0;
        500,0,15,270,0,0];
ustat2(:,2:6) = deg2rad(ustat2(:,2:6));   %degrees to radians                            
keps2(:,2:6) = deg2rad(keps2(:,2:6));    
ustat2(:,1)  = ustat2(:,1) + astroConstants(24);    %altitude to radius in [km]
keps2(:,1)   = keps2(:,1) + astroConstants(24); 
sit2 = 2;          %[-] 1 - Mars ground to Mars orbiter, %2 - Mars orbiter 
                   % to Mars orbiter, 3 - Mars to Earth (generic)
cas2 = 2;          %[-] case number
frq2 = 401.6e6;     %[Hz] carrier signal frequency
powt2 = 40;        %[W] ground user RF power emitted
BW2 = 1e6;         %[Hz] Bandwidth
point2 = [0 0];    %[tran recv] 1 - antenna steered, 0 - not steered
%Run the function and output graphs and totals
[UHF_Orb] = pass_over(1,sit2,cas2,frq2,powt2,elec,keps2,ustat2,t,dt,point2,...
                  UHF_G,UHF_G,'Orbiter to Orbiter UHF',BW2,1);

%% (PROOF) CASE 3 - MRO to Earth
ustat3 = [a_Mars,0,0,0,0,0];
keps3 = [a_Earth,0,0,0,0,0];     %[km,deg,deg,deg,deg,deg]
                                 %[a,e,i,OM,om,th,mu]
sit3 = 3;          %[-] 1 - Mars ground to Mars orbiter, %2 - Mars orbiter 
                   % to Mars orbiter, 3 - Mars to Earth (generic)
cas3 = 3;          %[-] case number
frq3 = 8440e6;     %[Hz] carrier signal frequency
powt3 = 100;        %[W] ground user RF power emitted
BW3 = 5e5;         %[Hz] Bandwidth
point3 = [1 1];    %[tran recv] 1 - antenna steered, 0 - not steered
%Run the function and output graphs and totals
[MRO_ME] = pass_over(0,sit3,cas3,frq3,powt3,dsn,keps3,ustat3,t,dt,point3,...
                  S_HGA_X,DSN_MGA_X,'Mars to Earth',BW3,1);

%% CASE 3 - Something to do with Phased array antennas ???