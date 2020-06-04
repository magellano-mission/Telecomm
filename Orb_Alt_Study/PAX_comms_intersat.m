
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


dt = 60; sols = 1; t = 0: dt : sols*88620; %[s] Mday=88620, Eday=86400
sit.cas = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
sit.pol = 0;
frq = 8490e6;    %[Hz] carrier signal frequency
powt = 30;        %[W] ground user RF power emitted

%Custom Antennas
custt.type = 'phased array';
custt.gain_peak = 27.4;
custt.HPBW = 2.7;
custt.dir = 1;
custt.tilt = [0 0];         %[deg] tilt from zenith in [azi ele] spherical directions
custt.plotting = 0;
custt.lims = [-60 60];
custr.type = 'phased array';
custr.gain_peak = 27.4;
custr.dir = -1;
custr.HPBW = 2.7;
custr.tilt = [0 0];         %[deg] tilt from zenith in [azi ele] spherical directions
custr.lims = [-60 60];
custr.plotting = 0;

hard = sys_hard(0,0,custt,custr,'sdst',290,[0 0]);

%% FUNCTION
tic; out = struct; res = NaN(length(orb_alts1),length(orb_alts2),8);

for i = 1:length(orb_alts1)
    for j = 1:length(orb_alts2)
        if orb_alts1(i) > orb_alts2(j)
            [out(i,j).all,tots] = pass_over(sit,frq,powt,hard,keps1(i,:),keps2(j,:),t,dt,0);
            %storing totals                        %[kJ]    [uJ]    [Gb]    vis     con 
            res(i,j,:) = [orb_alts1(i) orb_alts2(j) tots(1) tots(2) tots(3) tots(4) tots(5) sols];
        end
    end
end

S(1) = load('gong'); sound(S(1).y,S(1).Fs); toc

%% SAVE RESULTS
%save results and variables to output file for post-processing convenience
%save('Output Files\Intersat\PAX_comms_intersat_out','out','res','keps1','keps2',...
%     'dt','sols','sit','frq','powt','hard','t','orb_alts1','orb_alts2')

%% POST-PROCESSING
%Load information from file if required 
%load('Output Files\Intersat\PAX_comms_intersat_out')
plots = zeros(length(orb_alts1),length(orb_alts2),4);
plots(:,:,1) = res(:,:,5)./sols;             %[Gb/sol] Daily data transfer
plots(:,:,2) = res(:,:,5)./res(:,:,3)./sols;   %[Gb/kJ/sol] Transfer Efficiency
plots(:,:,3) = res(:,:,6)./sols;             %[hrs/sol] Visible Time
plots(:,:,4) = res(:,:,7)./sols;             %[hrs/sol] Downlink Time

%% STUDY PLOTTING

figure(1)
subplot(2,2,1)
colormap('winter')
levels = 100:50:500;
[C,h] = contour(orb_alts2,orb_alts1,plots(:,:,1),levels,'LabelSpacing',300);
clabel(C,h,'FontSize',10);
xlim([0 3000])
ylim([0 10000])
xlabel('RS Altitude [km]')
ylabel('ECS Altitude [km]')
title({'Orbiter Altitudes vs Daily Data Transfer [Gb]'})
grid on

subplot(2,2,2)
colormap('winter')
[C,h] = contour(orb_alts2,orb_alts1,plots(:,:,2)*1000,'LabelSpacing',300);
clabel(C,h,'FontSize',10);
xlim([0 3000])
ylim([0 10000])
xlabel('RS Altitude [km]')
ylabel('ECS Altitude [km]')
title('Orbiter Altitudes vs Data Transfer Efficiency [Mb/kJ/sol]')
grid on

%yyaxis left
%plot(res(:,1),plots(:,1))
%xlabel('ECS Altitude [km]')
%ylabel('Data Transfer [Gb/sol]')
%ylim([0 inf])
%title({'ECS Altitude vs Daily Data Transfer'})
%hold on
%yyaxis right
%plot(res(:,1),plots(:,2)*1000,'--')
%xlabel('ECS Altitude [km]')
%ylabel('Data Transfer Efficiency [Mb/kJ/sol]')
%ylim([0 inf])
%legend('Daily Data Transfer','Data Transfer Efficiency','Location','southwest')
%pbaspect([1 1 1])
%hold off

%subplot(1,2,2)
%plot(res(:,1),plots(:,3))
%xlabel('ECS Altitude [km]')
%ylabel('Visible Time  [hrs/sol]')
%title({'Orbiter Altitude vs Daily Time Windows'})
%hold on
%plot(res(:,1),plots(:,4),'--')
%ylim([0 inf])
%legend('Daily Visible Time','Daily Downlink Time','Location','southeast')
%pbaspect([1 1 1])
%hold off

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