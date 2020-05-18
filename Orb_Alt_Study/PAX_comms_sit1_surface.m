
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% INPUTS
orb_alts = 200:200:5000;              %[km] altitude range of interest
a = orb_alts + astroConstants(24);     %[km] semi-major axis range of interest
keps = zeros(length(a),6);              %[km & rads] 
keps(:,1) = 1*a';                       %[km & rads]
keps(:,3) = pi/9;

lats = -30:2.5:30;
ustat = zeros(length(lats),2);
ustat(:,1) = 1*lats';              %[deg lat, deg long] typical user position
ustat = deg2rad(ustat);            %degrees to radians

dt = 60; sols = 5; t = 0: dt : sols*88620; %[s] Mday=88620, Eday=86400
sit = 1;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
frq = 8490e6;    %[Hz] carrier signal frequency
powt = 15;        %[W] ground user RF power emitted

%Custom Antennas
custt.type = 'phased array';
custt.gain_peak = 18.2;
custt.HPBW = 17.1;
custt.plotting = 0;
custr.type = 'phased array';
custr.gain_peak = 26.8;
custr.HPBW = 6.4;
custr.plotting = 0;

hard = sys_hard(0,0,custt,custr,'sdst',290,[0 0]);

%% FUNCTION
tic; out = struct; res = zeros(length(lats),length(a),8);
for i = 1:length(lats)
    for j = 1:length(a)
    [out(i,j).all,tots] = pass_over(sit,frq,powt,hard,keps(j,:),ustat(i,:),t,dt,0);
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
plots = zeros(length(lats),length(out),6);
plots(:,:,1) = res(:,:,5)./sols;               %[Gb/sol] Daily data transfer
plots(:,:,2) = res(:,:,5)./res(:,:,3)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,:,3) = res(:,:,6)./sols;               %[hrs/sol] Visible Time
plots(:,:,4) = res(:,:,7)./sols;               %[hrs/sol] Downlink Time

%% STUDY PLOTTING
figure(1)
subplot(2,2,1)
surf(orb_alts,lats,plots(:,:,1))
ylabel('Gound User Latitude [km]')
xlabel('Orbiter Altitude [km]')
zlabel('Daily Data Transfer [Gb]')
zlim([0 inf])

subplot(2,2,2)
surf(orb_alts,lats,plots(:,:,2)*1000)
ylabel('Gound User Latitude [km]')
xlabel('Orbiter Altitude [km]')
zlabel('Data Transfer Efficiency [Mb/kJ/sol]')
zlim([0 inf])

subplot(2,2,3)
surf(orb_alts,lats,plots(:,:,3))
ylabel('Gound User Latitude [degrees]')
xlabel('Orbiter Altitude [km]')
zlabel('Daily Time Window [hrs]')
hold on
surf(orb_alts,lats,plots(:,:,4))
zlim([0 inf])
legend('Visible Time','Downlink Time','Location','southeast')
hold off

fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

%% BREAK
return

%% PARTICULAR CASE PLOTTING
%choose an altitude of interest to produce plots for
alt_int = 17000;     %[km] altitude
[log,row] = ismember(alt_int,orb_alts);

%choose a new time frame to produce the plots for
dt_plot = 10;
sols_plot = 1;
t_plot = 0:dt_plot:88620*sols_plot;

if log == 1
[~,~] = pass_over(sit,frq,powt,hard,keps(row,:),ustat,t_plot,dt_plot,1);
end