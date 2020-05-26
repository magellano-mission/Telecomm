
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% INPUTS
orb_alts = 4300:20:4500;              %[km] altitude range of interest
incs = -25:5:25;                        %[degrees] RS inclination range of interest
inc = deg2rad(incs);
a = orb_alts + astroConstants(24);     %[km] semi-major axis range of interest
keps = zeros(length(orb_alts),length(inc),6);              %[km & rads] 

for jj = 1:length(orb_alts)
    for ii = 1:length(incs)
        keps(jj,ii,1) = a(jj);
        keps(jj,ii,3) = inc(ii);
    end
end

ustat = [4900,0,deg2rad(-25),0,0,0];             %[km alt & degrees] typical user position
%ustat(1) = ustat(1) + astroConstants(24);
%ustat(2:6) = deg2rad(ustat(2:6));     %degrees to radians

dt = 60; sols = 5; t = 0: dt : sols*88620; %[s] Mday=88620, Eday=86400
sit = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
frq = 8490e6;    %[Hz] carrier signal frequency
freq = 'X band';
powt = 25;        %[W] ground user RF power emitted

%Custom Antennas
custt.type = 'phased array';
custt.gain_peak = 27.4;
custt.HPBW = 2.7;
custt.plotting = 0;
custr.type = 'phased array';
custr.gain_peak = 27.4;
custr.HPBW = 2.7;
custr.plotting = 0;

hard = sys_hard(0,0,custt,custr,'sdst',290,[0 0]);

%% FUNCTION
tic; out = struct; res = zeros(length(incs),length(a),8);
for i = 1:length(incs)
    for j = 1:length(orb_alts)
    [out(i,j).all,tots] = pass_over(sit,frq,powt,hard,keps(j,i,:),ustat,t,dt,0);
        %storing totals                    %[kJ]    [uJ]    [Gb]    vis     con 
        res(i,j,:) = [incs(i) orb_alts(j) tots(1) tots(2) tots(3) tots(4) tots(5) sols]; 
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
plots = zeros(length(incs),length(orb_alts),6);
plots(:,:,1) = res(:,:,5)./sols;               %[Gb/sol] Daily data transfer
plots(:,:,2) = res(:,:,5)./res(:,:,3)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,:,3) = res(:,:,6)./sols;               %[hrs/sol] Visible Time
plots(:,:,4) = res(:,:,7)./sols;               %[hrs/sol] Downlink Time

%% STUDY PLOTTING
figure(1)
subplot(2,2,1)
surf(orb_alts,incs,plots(:,:,1))
ylabel('RS Orbiter inclincation [degrees]')
xlabel('RS Orbiter Altitude [km]')
zlabel('Daily Data Transfer [Gb]')
zlim([0 inf])

subplot(2,2,2)
surf(orb_alts,incs,plots(:,:,2)*1000)
ylabel('RS Orbiter inclincation [degrees]')
xlabel('RS Orbiter Altitude [km]')
zlabel('Data Transfer Efficiency [Mb/kJ/sol]')
zlim([0 inf])

subplot(2,2,3)
surf(orb_alts,incs,plots(:,:,3),'FaceColor','b')
ylabel('RS Orbiter inclincation [degrees]')
xlabel('RS Orbiter Altitude [km]')
zlabel('Daily Time Window [hrs]')
hold on
surf(orb_alts,incs,plots(:,:,4),'FaceColor','r')
zlim([0 inf])
legend('Visible Time','Downlink Time','Location','southeast')
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