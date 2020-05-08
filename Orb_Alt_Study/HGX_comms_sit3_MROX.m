
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% INPUTS
orb_alt_m = 1.52*astroConstants(2);   %[km] Mars average orbital radius
keps_m = [orb_alt_m 0 0 0 0 0];       %[km & rads]

orb_alt_e = astroConstants(2);        %[km] Earth average orbital radius
keps_e = [orb_alt_e 0 0 0 0 0];       %[km & rads]

dt = 86400; days = 365*2.135; t = 0: dt : days*86400; %[s] Mday=88620, Eday=86400
sit = 3;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
frq = 8490e6;    %[Hz] carrier signal frequency
powt = 100;        %[W] use user RF power emitted

hard = sys_hard('OHGAX','DSN34X','dsn',290,[1 1]);
hard.cont.BW = 2e6;

%% FUNCTION
tic; out = struct; %res = zeros(length(keps_m),7);

[out.all,tots] = pass_over(sit,frq,powt,hard,keps_m,keps_e,t,dt,0);
    %storing totals          %[kJ]    [uJ]    [Gb]    vis     con 
    res = [keps_m(1) tots(1) tots(2) tots(3) tots(4) tots(5) days]; 


S(1) = load('gong'); sound(S(1).y,S(1).Fs); toc

%% SAVE RESULTS
%save results and variables to output file for post-processing convenience
%save('Output Files\PAX_comms_inter_system_out','out','res','ustat','keps',...
%     'dt','sols','sit','frq','powt','point','t')

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

%% BREAK
return

%% PARTICULAR CASE PLOTTING
%choose an altitude of interest to produce plots for
%alt_int = 17000;     %[km] altitude
%[log,row] = ismember(alt_int,orb_alts);

%choose a new time frame to produce the plots for
%dt_plot = 10;
%sols_plot = 1;
%t_plot = 0:dt_plot:88620*sols_plot;

%if log == 1
[~,~] = pass_over(sit,frq,powt,hard,keps_m,keps_e,t,dt,1);
%end