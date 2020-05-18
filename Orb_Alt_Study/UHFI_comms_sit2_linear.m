
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

%Environment parameters
T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space (REVISIT)
a_Earth = astroConstants(2);
a_Mars  = 1.52*a_Earth;

%% INPUTS - TYPICAL ORBITAL USER UHF STUDY
%Ground user is kept at fixed equatorial location with fixed hardware and 
%fixed power.Eclipses are not considered. The orbital altitude of the
%orbiter is varied
orb_alts = 200:100:10000;              %[km] altitude range of interest
a = orb_alts + astroConstants(24);     %[km] semi-major axis range of interest
keps = zeros(length(a),6);              %[km & rads] 
keps(:,1) = 1*a';                       %[km & rads]

ustat = [300,0,75,0,0,0];             %[km alt & degrees] typical user position
ustat(1) = ustat(1) + astroConstants(24);
ustat(2:6) = deg2rad(ustat(2:6));     %degrees to radians

dt = 15;                %[s] time step
sols = 5;               
t = 0: dt : sols*88620; %[s] Mars day = 88620, Earth day = 86400

sit = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter 
                  % to Mars orbiter, 3 - Mars to Earth (generic)
frq = 401.6e6;    %[Hz] carrier signal frequency
powt = 30;        %[W] user RF power emitted

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

%% FUNCTION - TYPICAL ORBITAL USER UHF STUDY
tic
out = struct;
res = zeros(length(a),7);
for i = 1:length(a)
    [out(i).all,tots] = pass_over(sit,frq,powt,hard,keps(i,:),ustat,t,dt,0);
    %storing totals          %[kJ]    [uJ]    [Gb]    vis     con 
    res(i,:) = [orb_alts(i) tots(1) tots(2) tots(3) tots(4) tots(5) sols];
end
S(1) = load('chirp');
sound(S(1).y,S(1).Fs)
toc

%% SAVE RESULTS
%save results and variables to output file for post-processing convenience
%save('Output Files\UHF_comms_trends_sit2_out','out','res','ustat','keps',...
%     'dt','sols','sit','frq','powt','point','t')

%% POST-PROCESSING
%Load information from file if required 
%load('Output Files\UHF_comms_trends_sit2_out')
plots = zeros(length(out),6);
plots(:,1) = res(:,4)./sols;             %[Gb/sol] Daily data transfer
plots(:,2) = res(:,4)./res(:,2)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,3) = res(:,5)./sols;             %[hrs/sol] Visible Time
plots(:,4) = res(:,6)./sols;             %[hrs/sol] Downlink Time


%% STUDY PLOTTING
figure(1)
subplot(1,2,1)
yyaxis left
plot(res(:,1),plots(:,1))
xlabel('Altitude [km]')
ylabel('Data Transfer [Gb/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Daily Data Transfer'})
hold on
yyaxis right
plot(res(:,1),plots(:,2)*1000,'--')
xlabel('Altitude [km]')
ylabel('Data Transfer Efficiency [Mb/kJ/sol]')
ylim([0 inf])
legend('Daily Data Transfer','Data Transfer Efficiency','Location','southwest')
pbaspect([1 1 1])
hold off

subplot(1,2,2)
plot(res(:,1),plots(:,3))
xlabel('Altitude [km]')
ylabel('Visible Time  [hrs/sol]')
title({'Orbiter Altitude vs Daily Time Windows'})
hold on
plot(res(:,1),plots(:,4),'--')
ylim([0 inf])
legend('Daily Visible Time','Daily Downlink Time','Location','southeast')
pbaspect([1 1 1])
hold off

fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

%% BREAK
return

%% PARTICULAR CASE PLOTTING 
%choose an altitude of interest to produce plots for
alt_int = 3000;     %[km] altitude
[log,row] = ismember(alt_int,orb_alts);

%choose a new time frame to produce the plots for
dt_plot = 10;
sols_plot = 1;
t_plot = 0:dt_plot:88620*sols_plot;

if log == 1
[~,~] = pass_over(sit,frq,powt,hard,keps(row,:),ustat,t_plot,dt_plot,1);
end