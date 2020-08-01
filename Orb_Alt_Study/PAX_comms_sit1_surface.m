
%% SET UP
clearvars; clc; %close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% INPUTS
orb_alts = 200:200:3000;              %[km] altitude range of interest
a = orb_alts + astroConstants(24);     %[km] semi-major axis range of interest
keps = zeros(length(a),6);              %[km & rads] 
keps(:,1) = 1*a';                       %[km & rads]
inc = 0;
keps(:,3) = deg2rad(inc);   % 15 degree inclincation

lats = 0:3:30;
ustat = zeros(length(lats),2);
ustat(:,1) = 1*lats';              %[deg lat, deg long] typical user position
ustat = deg2rad(ustat);            %degrees to radians

dt = 30; sols = 5; t = 0: dt : sols*88620; %[s] Mday=88620, Eday=86400
sit.cas = 1;   %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
sit.ops = 1;   %[-] 0 - nominal conditions, 1 - worst case conditions 
sit.pol = 0;
frq = 8490e6;    %[Hz] carrier signal frequency
freq = 'X band';
powt = 10;        %[W] ground user RF power emitted

%Custom Antennas
custt.type = 'phased array';
custt.gain_peak = 21.3;
custt.dir = 1;
custt.HPBW = 2.7;
custt.tilt = [0 0];         %[deg] tilt from zenith in [azi ele] spherical directions
custt.lims = [-60 60];
custt.plotting = 0;
custr.type = 'phased array';
custr.gain_peak = 24.3;
custr.dir = -1;
custr.HPBW = 2.7;           %[deg] tilt from nadir in [azi ele] spherical directions
custr.tilt = [0 -40];
custr.lims = [-60 60];
custr.plotting = 0;

hard = sys_hard(0,0,custt,custr,'sdst',290,[0 0],sit);

%% FUNCTION
tic; out = struct; res = zeros(length(lats),length(a),8);
for i = 1:length(lats)
    for j = 1:length(a)
    [out(i,j).all,~,tots] = pass_over(sit,frq,powt,hard,keps(j,:),ustat(i,:),t,dt,0);
        %storing totals                    %[kJ]    [uJ]    [Gb]    vis     con 
        res(i,j,:) = [ustat(i) orb_alts(j) tots(1) tots(2) tots(3) tots(4) tots(5) sols]; 
    end
end
S(1) = load('chirp'); sound(S(1).y,S(1).Fs); toc

%% SAVE RESULTS
%save results and variables to output file for post-processing convenience
%save('Output Files\PAX_comms_inter_system_out','out','res','ustat','keps',...
%     'dt','sols','sit','frq','powt','point','t')

%% POST-PROCESSING
%Load information from file if required 
%load('Output Files\???')
plots = zeros(length(lats),length(orb_alts),6);
plots(:,:,1) = res(:,:,5)./sols;               %[Gb/sol] Daily data transfer
plots(:,:,2) = res(:,:,5)./res(:,:,3)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,:,3) = res(:,:,6)./sols;               %[hrs/sol] Visible Time
plots(:,:,4) = res(:,:,7)./sols;               %[hrs/sol] Downlink Time

%% STUDY PLOTTING
%% Adding to path
addpath(genpath(fileparts(pwd)))

%% Figure Initialization    
load('MagellanoColorMap.mat');
colormap(MagellanoColorMap)

figure(1)
subplot(2,2,1)
levels =  [10 20 30 40 50 60 80 100];
[C,h] = contour(orb_alts,lats,plots(:,:,1),levels);
clabel(C,h,'manual','BackgroundColor','white','Margin',0.5,'FontSize',10)
h.LineWidth = 2.5;
grid on
ylabel('Gound User Latitude [degrees]')
ylim([0 30])
xlabel('Orbiter Altitude [km]')
title('Daily Data Transfer [Gb]')
zlim([0 inf])

subplot(2,2,2)
levels =  [50 100 150 200 300];
[C,h] = contour(orb_alts,lats,plots(:,:,2)*1000,levels);
clabel(C,h,'manual','BackgroundColor','white','Margin',0.5,'FontSize',10)
h.LineWidth = 3;
grid on
ylabel('Gound User Latitude [degrees]')
ylim([0 30])
xlabel('Orbiter Altitude [km]')
title('Data Transfer Efficiency [Mb/kJ/sol]')
zlim([0 inf])

subplot(2,2,3)
levels = 0:1:24;
[C,h] = contour(orb_alts,lats,plots(:,:,3),levels);
clabel(C,h,'manual','BackgroundColor','white','Margin',0.5,'FontSize',10)
h.LineWidth = 3;
grid on
ylabel('Gound User Latitude [degrees]')
ylim([0 30])
xlabel('Orbiter Altitude [km]')
title('Daily Visibility to RS [hrs]')

subplot(2,2,4)
levels = 0:1:24;
[C,h] = contour(orb_alts,lats,plots(:,:,4),levels);
clabel(C,h,'manual','BackgroundColor','white','Margin',0.5,'FontSize',10)
h.LineWidth = 3;
grid on
ylabel('Gound User Latitude [degrees]')
ylim([0 30])
xlabel('Orbiter Altitude [km]')
title('Daily Downlink to RS Availability [hrs]')

sgtitle({powt+"W (RF) User with "+custt.gain_peak+"dBi peak gain "+freq+" antenna", ...
        inc+" degree inclination orbiter with "+custr.gain_peak+"dBi gain "+freq+" antenna"}) 

fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

%hold on
%plot(lats,plots(:,1,1))
%xlabel('Surface User Latitude [degrees]')
%ylabel('Data Return [Gb]')
%grid on


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