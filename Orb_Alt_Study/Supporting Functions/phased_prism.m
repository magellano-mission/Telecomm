function [prism] = phased_prism(sides,peak_gain,plotting)

% This function takes the number of sides of the ECS phased array prism and
% the peak gain and calculates the antenna gain at all places on a sphere
% around the antenna. The phased array steering losses are calculated
% assuming the n factor of 2.4 from the MESSENGER phased arrays. Plots of
% the 2D and 3D gain patterns are output if requested. The elev = 0 gain
% pattern is output as a matrix for later interpolation.

% --- TO DO ---
% Investigate outputting the full 3D gain surface for better accuracy of
% gain estimation for out of plane communication


%% INPUTS
% sides -       [-] number of sides on the prism, must be 2 or more
% peak_gain -   [dBi] peak gain of each phased array
% plotting -    [-] 0 for no plot, 1 for plot

%retrieve the antenna gain curve
type = 'phased array';
HPBW = NaN;
[curve] = basic_antenna(type,peak_gain,HPBW,0);

%% SET UP
%Create the geometry of the prism
pris_azi = linspace(0,2*pi,sides+1);    %[rad] surface azimuth normals (sides+1)
pris_ele = zeros(1,length(pris_azi));  %[rad] surface elevation = 0       
pris_rad = 1;                           %[-] unit radius

%Create array of surface Cartesian normal vectors
[pris(:,1),pris(:,2),pris(:,3)] = sph2cart(pris_azi,pris_ele,pris_rad);
pris = pris(1:end-1,:);

%Create points around the global perimeter to evaluate the gain
azi = linspace(0,2*pi,100);      %vector of azimuths
ele_deg = -90:5:90;
ele = deg2rad(ele_deg);  %vector of elevations

%% LOOP
%initialise blank matrices to store gains
gains3 = NaN(length(azi),length(ele),sides+1);
gains2 = NaN(length(azi),length(ele),sides+1);

for i = 1:length(ele)       %loop over all elevations
    for j =1:length(azi)    %loop over all azimuths
        % cartesian normalised position vector
        [pnt(1),pnt(2),pnt(3)] = sph2cart(azi(j),ele(i),1);
        for k = 1:sides     %loop over all array panels
            % Angle between position vector and array surface normal
            z = atan(norm(cross(pris(k,:),pnt))/dot(pris(k,:),pnt));
            % get gain if angle is within functional range
            if abs(z) <= pi/3 && z >= 0
                gains3(j,i,k) = curve.bor(z);
                gains2(j,i,k) = gains3(j,i,k);
            end
        end
        
        % use the highest gain offered by any of the panels
        gains3(j,i,end) = nanmax(gains3(j,i,1:(end-1)));
        gains2(j,i,end) = gains3(j,i,end);
        
        % eliminate values which cause problems for the plotting functions
        if abs(gains3(j,i,end)) == Inf || isnan(gains3(j,i,end))
            gains3(j,i,end) = 0;
            gains2(j,i,end) = NaN;
        end
        if gains3(j,i,end) < 2
            gains3(j,i,end) = 0;
            gains2(j,i,end) = NaN;
        end
    end
end

%output the 2D array of ele=0 azimuth vs gain for ECS to RS communications
temp = gains2(:,:,end);
prism(:,1) = azi';
prism(:,2) = temp(:,19);

%% PLOTS
if plotting == 1
    
% plot the 3D gain pattern
figure(8)
hplot = patternCustom(gains3(:,:,end),90-rad2deg(ele),rad2deg(azi));
caxis([peak_gain-0.3*peak_gain, peak_gain])
bar = colorbar;
ylabel(bar,'Antenna Gain [dBi]')
set(hplot,'LineStyle','-');
axis equal
title('Phased Array Prism 3D Gain Pattern')


%plot the 2D gain curves, paramaterised by elevation
figure(9)
rot = 0;
for m = 1:length(ele_deg)
    r = rem(ele_deg(m),15);
    if r == 0 && ele_deg(m) >= 0 && ele_deg(m) <= 50
        rot = rot + 1;
        polarplot(azi,temp(:,m),'DisplayName',['elev.= +/- ' num2str(ele_deg(m)) ' degrees'],'LineWidth',0.8)
        hold on
    end
end
legend('show','Position', [0.5 0.5 0 0])
title('Prism Gain [dBi] vs Azimuth and Elevation [degrees]')
rticks(0:3:30)
hold off
end
