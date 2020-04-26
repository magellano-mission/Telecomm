%Later this script will get turned into a function called by
%system_characteristics
clearvars
clc
%% Earth - Mars ECS
%Shared properties
frq = 34.5e9;   %[Hz] Signal nominal frequency, 8 GHz for MRO
BW = 1.75e8;     %[Hz] nominal bandwidth
dist = 400e9;  %[m] Distance between transmitter and receiver
%design margin

%Transmitter properties
pow_t = 2e5;   %[W] Signal power sent to antenna, 20kW for X band DSN
dia_t = 34;    %[m] Diameter of transmitting dish
e_a = 0.94;   %[-] Aperture efficiency of transmitting antenna
atm_t = 0.9;    %[-] Transmission coefficient through atmosphere (EARTH)
%rain fade

%Receiver properties
dia_r = 4.5;    %[m] Diameter of receiving dish
k_r = 0.8;      %[-] Aperture efficiency of receiving dish

Ear_ECS.name = 'Earth to Mars ECS Uplink Budget';
[Ear_ECS.budget] = zfun_link_budget(frq,BW,dist,pow_t,dia_t,e_a,atm_t, ...
    dia_r,k_r);

disp(Ear_ECS.name)
disp(Ear_ECS.budget)

%% Mars ECS - Mars RS
% Call the link budget function

%% Mars RS - Mars RS
% Call the link budget function

%% Mars RS - Mars Surface User
% Call the link budget function

%% Mars RS - Mars Orbital User
% Call the link budget function

%Later add off-nominal cases when Earth or ECS talk directly to Users in
%orbit or on surface. Later can consider if we want to add additional links
%in the chain back to Earth or around Mars



