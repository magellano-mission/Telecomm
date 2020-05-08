
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% INPUTS
orb_alts1 = 5000:500:17000;              %[km] altitude range of interest
a1 = orb_alts1 + astroConstants(24);     %[km] semi-major axis range of interest
keps1 = zeros(length(a1),6);              %[km & rads] 
keps1(:,1) = 1*a1';                       %[km & rads]
keps1(:,3) = 0;

orb_alts2 = 200:100:5000;              %[km] altitude range of interest
a2 = orb_alts2 + astroConstants(24);     %[km] semi-major axis range of interest
keps2 = zeros(length(a2),6);              %[km & rads] 
keps2(:,1) = 1*a2';                       %[km & rads]
keps2(:,3) = 0;


dt = 30; sols = 1; t = 0: dt : sols*88620; %[s] Mday=88620, Eday=86400
sit = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
frq = 8490e6;    %[Hz] carrier signal frequency
powt = 40;        %[W] ground user RF power emitted

hard = sys_hard('PAX','PAX','sdst',290,[0 0]);

%% FUNCTION
tic; out = struct; res = zeros(length(a1),length(a2),8);

for i = 1:length(a1)
    for j = 1:length(a2)
        [out(i,j).all,tots] = pass_over(sit,frq,powt,hard,keps1(i,:),keps2(j,:),t,dt,0);
        %storing totals                      %[kJ]    [uJ]    [Gb]    vis     con 
        res(i,j,:) = [orb_alts1(i) orb_alts2(j) tots(1) tots(2) tots(3) tots(4) tots(5) sols];
    end
end

S(1) = load('gong'); sound(S(1).y,S(1).Fs); toc

%% SAVE RESULTS
%save results and variables to output file for post-processing convenience
save('Output Files\PAX_comms_intersat_out','out','res','keps1','keps2',...
     'dt','sols','sit','frq','powt','hard','t')

%% POST-PROCESSING
%Load information from file if required 
%load('Output Files\???')
plots = zeros(i,j,4);
plots(:,:,1) = res(:,:,5)./sols;             %[Gb/sol] Daily data transfer
plots(:,:,2) = res(:,:,5)./res(:,:,3)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,:,3) = res(:,:,6)./sols;             %[hrs/sol] Visible Time
plots(:,:,4) = res(:,:,7)./sols;             %[hrs/sol] Downlink Time

%% STUDY PLOTTING
figure(1)
subplot(2,2,1)
surf(orb_alts2,orb_alts1,plots(:,:,1)*1000)
xlabel('Orbiter 1 Altitude [km]')
ylabel('Orbiter 2 Altitude [km]')
zlabel('Data Transfer [Mb/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Daily Data Transfer'})

subplot(2,2,2)
surf(orb_alts2,orb_alts1,plots(:,:,2)*1000)
xlabel('Orbiter 1 Altitude [km]')
xlabel('Orbiter 2 Altitude [km]')
zlabel('Data Transfer Efficiency [Mb/kJ/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Data Transfer Efficiency'})

subplot(2,2,3)
surf(orb_alts2,orb_alts1,plots(:,:,3))
xlabel('Orbiter 1 Altitude [km]')
xlabel('Orbiter 2 Altitude [km]')
zlabel('Visible Time  [hrs/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Daily Visible Time'})

subplot(2,2,4)
surf(orb_alts2,orb_alts1,plots(:,:,4))
xlabel('Orbiter 1 Altitude [km]')
xlabel('Orbiter 2 Altitude [km]')
zlabel('Downlink Time  [hrs/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Daily Downlink Time'})

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
[~,~] = pass_over(sit,frq,powt,hard,keps1(row,:),ustat,t_plot,dt_plot,1);
end