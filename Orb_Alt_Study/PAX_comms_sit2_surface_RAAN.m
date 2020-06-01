
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% INPUTS
RAAN = linspace(0,2*pi,100);
keps = zeros(length(RAAN),6);
keps(:,1) = 7400;
keps(:,3) = deg2rad(-25);
keps(:,4) = 1.*RAAN';

ustat = [4900,0,deg2rad(-25),0,0,0];             %[km alt & degrees] typical user position

dt = 60; sols = 5; t = 0: dt : sols*88620; %[s] Mday=88620, Eday=86400
sit = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
frq = 8490e6;    %[Hz] carrier signal frequency
freq = 'X band';
powt = 30;        %[W] RF power emitted

%Custom Antennas
custt.type = 'phased array';
custt.gain_peak = 27.4;
custt.HPBW = 2.7;
custt.tilt = [0 0];         %[deg] tilt from zenith in [ele azi] spherical directions
custt.plotting = 0;
custr.type = 'phased array';
custr.gain_peak = 27.4;
custr.HPBW = 2.7;
custr.tilt = [0 0];         %[deg] tilt from zenith in [ele azi] spherical directions
custr.plotting = 0;

hard = sys_hard(0,0,custt,custr,'sdst',290,[0 0]);

%% FUNCTION
tic; out = struct; res = zeros(length(RAAN),6);

    for j = 1:length(RAAN)
    [out(j).all,tots] = pass_over(sit,frq,powt,hard,keps(j,:),ustat,t,dt,0);
        %storing totals    %[kJ]    [uJ]    [Gb]    vis     con 
        res(j,:) = [tots(1) tots(2) tots(3) tots(4) tots(5) sols]; 
    end

S(1) = load('chirp'); sound(S(1).y,S(1).Fs); toc

%% SAVE RESULTS
%save results and variables to output file for post-processing convenience
%save('Output Files\PAX_comms_inter_system_out','out','res','ustat','keps',...
%     'dt','sols','sit','frq','powt','point','t')

%% POST-PROCESSING
%Load information from file if required 
%load('Output Files\???')
plots = zeros(length(RAAN),6);
plots(:,1) = res(:,3)./sols;               %[Gb/sol] Daily data transfer
plots(:,2) = res(:,3)./res(:,1)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,3) = res(:,4)./sols;               %[hrs/sol] Visible Time
plots(:,4) = res(:,5)./sols;               %[hrs/sol] Downlink Time

%% STUDY PLOTTING
figure(1)
subplot(2,2,1)
plot(rad2deg(RAAN),plots(:,1))
ylabel('Daily Data Transfer [Gb]')
xlabel('RAAN Shift [deg]')
xlim([0 360])
ylim([0 inf])
title('RS to ECS Daily Data Transfer')
grid on

subplot(2,2,2)
plot(rad2deg(RAAN),plots(:,2)*1000)
xlabel('RAAN Shift [deg]')
ylabel('Data Transfer Efficiency [Mb/kJ/sol]')
xlim([0 360])
ylim([0 inf])
title('RS to ECS Daily Data Transfer Efficiency')
grid on

subplot(2,2,3)
plot(rad2deg(RAAN),plots(:,3))
xlabel('RAAN Shift [deg]')
ylabel('Daily Time Window [hrs]')
hold on
plot(rad2deg(RAAN),plots(:,4))
legend('Visible Time','Downlink Time','Location','southeast')
title('RS to ECS Daily Time Windows')
xlim([0 360])
grid on
hold off

sgtitle({powt+"W (RF), 75 degree inclination orbital user with "+custt.gain_peak+"dBi peak gain "+freq+" antenna", ...
        "RS orbiter with "+custr.gain_peak+"dBi gain "+freq+" antenna"}) 

fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

%% BREAK
return

%% PARTICULAR CASE PLOTTING
%choose an altitude of interest to produce plots for
alt_int = 3000;     %[km] altitude
[~,row1] = ismember(alt_int,orb_alts);
keps_int = keps(row1,:).* ones(4,1);
keps_int(:,6) = [0;pi/2;pi;3*pi/2];

lats_int = 30;
[~,row2] = ismember(deg2rad(lats_int),ustat);

%choose a new time frame to produce the plots for
dt_plot = 10;
sols_plot = 1;
t_plot = 0:dt_plot:88620*sols_plot;

[~,~] = pass_over(sit,frq,powt,hard,keps_int,ustat(row2,:),t_plot,dt_plot,1);