function [] = sys_anal(orbstat,ustat,time,sit,hard)
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% FUNCTION
tic; out = struct; res = zeros(length(a),7);

for i = 1:length(a)
[out(i).all,tots] = pass_over(sit,frq,powt,hard,keps(i,:),ustat,t,dt,0);
    %storing totals          %[kJ]    [uJ]    [Gb]    vis     con 
    res(i,:) = [orb_alts(i) tots(1) tots(2) tots(3) tots(4) tots(5) sols]; 
end

S(1) = load('gong'); sound(S(1).y,S(1).Fs); toc

%% POST-PROCESSING
%Load information from file if required 
%load('Output Files\???')
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
ylim([0 inf])
title({'Orbiter Altitude vs Daily Data Transfer'})

subplot(2,2,2)
plot(res(:,1),plots(:,2)*1000)
xlabel('Altitude [km]')
ylabel('Data Transfer Efficiency [Mb/kJ/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Data Transfer Efficiency'})

subplot(2,2,3)
plot(res(:,1),plots(:,3))
xlabel('Altitude [km]')
ylabel('Visible Time  [hrs/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Daily Visible Time'})

subplot(2,2,4)
plot(res(:,1),plots(:,4))
xlabel('Altitude [km]')
ylabel('Downlink Time  [hrs/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Daily Downlink Time'})

fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

%% PARTICULAR CASE PLOTTING
return
%choose an altitude of interest to produce plots for
alt_int = 10000;     %[km] altitude
[log,row] = ismember(alt_int,altitudes);

%choose a new time frame to produce the plots for
dt_plot = 10;
sols_plot = 1;
t_plot = 0:dt_plot:88620*sols_plot;

if log == 1
[a,~] = pass_over(sit,2,frq,powt,sdst,keps(row,:),ustat,t_plot,dt_plot,point,1);
end