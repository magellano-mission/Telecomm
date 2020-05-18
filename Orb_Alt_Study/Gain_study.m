%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')

%% NOTES
% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% Environment parameters
T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space (REVISIT)
dt = 30;                  %[s] time step
t = 0: dt : 1*88620;      %[s] Mars day = 88620, Earth day = 86400
a_Earth = astroConstants(2);
a_Mars  = 1.52*a_Earth;


%% (PROOF) CASE 1 - MRO/Curiosity UHF at 401.6 MHz

%ustat1 = [0,0];
%ustat1 = deg2rad(ustat1);    %degrees to radians

ustat1 = [300,0,75,0,0,0];             %[km alt & degrees] typical user position
ustat1(1) = ustat1(1) + astroConstants(24);
ustat1(2:6) = deg2rad(ustat1(2:6));     %degrees to radians

keps1 = [6390,0,0,0,0,0;
         6390,0,0,0,0,pi/2;
         6390,0,0,0,0,pi;
         6390 0,0,0,0,3*pi/2];    %[altitude km,deg,deg,deg,deg,deg]
                             %[a,e,i,OM,om,th]

%keps1(:,2:6) = deg2rad(keps1(:,2:6));           %degrees to radians
%keps1(:,1)   = keps1(:,1) + astroConstants(24); %altitude to radius in [km]
sit1 = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter 
                   % to Mars orbiter, 3 - Mars to Earth (generic)
frq1 = 401.6e6;    %[Hz] carrier signal frequency
powt1 = 30;        %[W] ground user RF power emitted

%Custom Antennas
custt.type = 'helical';
custt.gain_peak = 8;
custt.HPBW = 60;
custt.plotting = 0;
custr.type = 'helical';
custr.gain_peak = 8;
custr.HPBW = 60;
custr.plotting = 0;

hard = sys_hard(0,0,custt,custr,'elec',290,[0 0]);

%Run the function and output graphs and totals
tic
[Cur_UHF] = pass_over(sit1,frq1,powt1,hard,keps1,ustat1,t,dt,1);
toc

%% (PROOF) CASE 3 - MRO to Earth
%ustat3 = [a_Mars,0,0,0,0,0];
%keps3 = [a_Earth,0,0,0,0,0];     %[km,deg,deg,deg,deg,deg]
                                 %[a,e,i,OM,om,th,mu]
%sit3 = 3;          %[-] 1 - Mars ground to Mars orbiter, %2 - Mars orbiter 
                   % to Mars orbiter, 3 - Mars to Earth (generic)
%cas3 = 3;          %[-] case number
%frq3 = 8440e6;     %[Hz] carrier signal frequency
%powt3 = 100;        %[W] ground user RF power emitted
%BW3 = 5e5;         %[Hz] Bandwidth
%point3 = [1 1];    %[tran recv] 1 - antenna steered, 0 - not steered
%Run the function and output graphs and totals
%[MRO_ME] = pass_over(0,sit3,cas3,frq3,powt3,dsn,keps3,ustat3,t,dt,point3,...
                  %S_HGA_X,DSN_MGA_X,'Mars to Earth',BW3,1);

%% CASE 3 - Something to do with Phased array antennas ???