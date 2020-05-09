
%% SET UP
clearvars; clc; close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')
addpath('Output Files')

% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% INPUTS
%Receiving orbiter in higher orbit
orb_alts1 = 200:500:17000;              %[km] altitude range of interest
a1 = orb_alts1 + astroConstants(24);     %[km] semi-major axis range of interest
keps1 = zeros(length(a1),6);              %[km & rads] 
keps1(:,1) = 1*a1';                       %[km & rads]
keps1(:,3) = 0;

%Transmitting orbiter in lower orbit
orb_alts2 = 200:200:10000;              %[km] altitude range of interest
a2 = orb_alts2 + astroConstants(24);     %[km] semi-major axis range of interest
keps2 = zeros(length(a2),6);              %[km & rads] 
keps2(:,1) = 1*a2';                       %[km & rads]
keps2(:,3) = 0;


dt = 30; sols = 3; t = 0: dt : sols*88620; %[s] Mday=88620, Eday=86400
sit = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
frq = 8490e6;    %[Hz] carrier signal frequency
powt = 40;        %[W] ground user RF power emitted

hard = sys_hard('UHFB','UHFB','elec',290,[0 0]);

%% FUNCTION
tic; out = struct; res = NaN(length(orb_alts1),length(orb_alts2),8);

for i = 1:length(orb_alts1)
    for j = 1:length(orb_alts2)
        if orb_alts1(i) > orb_alts2(j)
            [out(i,j).all,tots] = pass_over(sit,frq,powt,hard,keps1(i,:),keps2(j,:),t,dt,0);
            %storing totals                      %[kJ]    [uJ]    [Gb]    vis     con 
            res(i,j,:) = [orb_alts1(i) orb_alts2(j) tots(1) tots(2) tots(3) tots(4) tots(5) sols];
        end
    end
end

S(1) = load('gong'); sound(S(1).y,S(1).Fs); toc

%% SAVE RESULTS
%save results and variables to output file for post-processing convenience
%save('Output Files\Intersat\UHFB_comms_intersat_out','out','res','keps1','keps2',...
%     'dt','sols','sit','frq','powt','hard','t','orb_alts1','orb_alts2')

%% POST-PROCESSING
%Load information from file if required 
%load('Output Files\Intersat\UHFB_comms_intersat_out')
plots = zeros(length(orb_alts1),length(orb_alts2),4);
plots(:,:,1) = res(:,:,5)./sols;             %[Gb/sol] Daily data transfer
plots(:,:,2) = res(:,:,5)./res(:,:,3)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,:,3) = res(:,:,6)./sols;             %[hrs/sol] Visible Time
plots(:,:,4) = res(:,:,7)./sols;             %[hrs/sol] Downlink Time

%% STUDY PLOTTING
figure(1)
subplot(2,2,1)
surf(orb_alts2,orb_alts1,plots(:,:,1)*1000)
xlabel('RS Altitude [km]')
ylabel('ECS Altitude [km]')
zlabel('Data Transfer [Mb/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Daily Data Transfer'})

subplot(2,2,2)
surf(orb_alts2,orb_alts1,plots(:,:,2)*1000)
xlabel('RS Altitude [km]')
ylabel('ECS Altitude [km]')
zlabel('Data Transfer Efficiency [Mb/kJ/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Data Transfer Efficiency'})

subplot(2,2,3)
surf(orb_alts2,orb_alts1,plots(:,:,3))
xlabel('RS Altitude [km]')
ylabel('ECS Altitude [km]')
zlabel('Visible Time  [hrs/sol]')
ylim([0 inf])
title({'Orbiter Altitude vs Daily Visible Time'})

subplot(2,2,4)
surf(orb_alts2,orb_alts1,plots(:,:,4))
xlabel('RS Altitude [km]')
ylabel('ECS Altitude [km]')
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
alt_int1 = 10200;     %[km] altitude
[log1,row1] = ismember(alt_int1,orb_alts1);

alt_int2 = 1500;     %[km] altitude
[log2,row2] = ismember(alt_int2,orb_alts2);

%choose a new time frame to produce the plots for
dt_plot = 10;
sols_plot = 1;
t_plot = 0:dt_plot:88620*sols_plot;

if log1 == 1 && log2 == 1
[~,~] = pass_over(sit,frq,powt,hard,keps1(row1,:),keps2(row2,:),t_plot,dt_plot,1);
end