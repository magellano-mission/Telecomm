%% SET UP
clearvars; clc; %close all;
format compact
addpath('Antenna Gain Curves')
addpath('Supporting Functions')

%% NOTES
% Have included some basic dB margins in zfun_link_rate. They need to be
% revisited later to check accuracy and suitability.

%% Environment parameters
T_amb = 290;       %[K] Assumed ambient temperature (STANDARD - REVISIT)
T_ant = 20;        %[K] Assumed antenna temperature in space (REVISIT)
dt = 60;                  %[s] time step
t = 0: dt : 0.5*88620;      %[s] Mars day = 88620, Earth day = 86400
a_Earth = astroConstants(2);
a_Mars  = 1.52*a_Earth;


%% INPUTS
rstat = [4900,0,0,0,0,0];            
tstat = [4900,0,deg2rad(90),0,0,0];
%tstat = [deg2rad(0) 0];

sit1.cas = 2;          %[-] 1 - Mars ground to Mars orbiter, 2 - Mars orbiter to Mars orbiter, 3 - Mars to Earth (generic)
sit1.ops = 1;           %[-] 0 - nominal conditions, 1 - worst case conditions

frq1 = 8400e6;    %[Hz] carrier signal frequency
powt1 = 30;        %[W] ground user RF power emitted

%Custom Antennas
custt.type = 'phased array';
custt.gain_peak = 24.3;     %[dBi]
custt.HPBW = 2.7;
custt.dir = -1;              %[-] 1 zenith pointing, -1 nadir pointing
custt.tilt = [0 -45];         %[deg] tilt in [azim elev] spherical directions
custt.lims = [-60 60];      %[deg] steering limits
custt.plotting = 0;
custr.type = 'phased array';
custr.gain_peak = 24.3;     %[-] 1 zenith pointing, -1 nadir pointing
custr.HPBW = 2.7;
custr.dir = -1;             %[-] 1 zenith pointing, -1 nadir pointing
custr.tilt = [0 45];         %[deg] tilt in [azim elev] spherical directions
custr.lims = [-60 60];      %[deg] steering limits
custr.plotting = 0;

hard = sys_hard(0,0,custt,custr,'sdst',290,[0 0],sit1);

%% FUNCTION
tic
[recv,trans,tots] = pass_over(sit1,frq1,powt1,hard,rstat,tstat,t,dt,0);
toc

myVideo = VideoWriter('newfile.avi','Motion JPEG AVI');
myVideo.FrameRate = 10;
open(myVideo);

figure (1)
hold on
[x,y,z] = sphere(25);
surf(x*3390*1000,y*3390*1000,z*3390*1000,'facecolor',[0.9290, 0.6940, 0.1250]);
xlim([-4900000 4900000]);
ylim([-4900000 4900000]);
zlim([-4900000 4900000]);
set(gca,'visible','off')
axis equal


 for i = 1:length(recv.st.pos)
     pol(i) = plot3(recv.st.pos(i,1),recv.st.pos(i,2),recv.st.pos(i,3),'bo','MarkerSize',2);
     eq(i) = plot3(trans.st.pos(i,1),trans.st.pos(i,2),trans.st.pos(i,3),'ro','MarkerSize',2);
     y = plot3([recv.st.pos(i,1) trans.st.pos(i,1)], [recv.st.pos(i,2) trans.st.pos(i,2)], [recv.st.pos(i,3) trans.st.pos(i,3)],'k','LineWidth',1);
     view([12e6 3e6 3e6])
     frame = getframe(gcf);
     writeVideo(myVideo,frame);
     %pause(0.01)  
     if i < length(recv.st.pos)
        delete(y)
     end
     if i > 100
         delete(pol(1:i-100))
         delete(eq(1:i-100))
     end

 end
 %plot3(recv.st.pos(i,1),recv.st.pos(i,2),recv.st.pos(i,3),'ob','MarkerSize',5,'LineWidth',5)
 %plot3(trans.st.pos(i,1),trans.st.pos(i,2),trans.st.pos(i,3),'or','MarkerSize',5,'LineWidth',5)
 

