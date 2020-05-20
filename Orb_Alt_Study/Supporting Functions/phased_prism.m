function [prism] = phased_prism(sides,peak_gain,plotting)

% This function creates a 2D curve and table of the gain vs orientation for
% an ECS receiver with a prism of flat-panel phased arrays


%% INPUTS
% sides - [-] number of sides on the prism, must be 2 or more
% peak_gain - [dBi] peak gain of each phased array
% plotting - [-] 0 for no plot, 1 for plot

%Create the geometry of the prism
angles = linspace(0,2*pi,sides+1);      %surface normals

%Create points around the global perimeter to evaluate the gain
possy = linspace(0,2*pi,1000);
gains = zeros(1,length(possy));

%retrieve the antenna gain curve
type = 'phased array';
HPBW = NaN;
[curve] = basic_antenna(type,peak_gain,HPBW,0);

for i = 1:length(possy)     %iterate for each point on the circumference
    for j = 1:length(angles) %iterate through all the sides
        if possy(i) >= angles(j) - 0.33*pi ...
        && possy(i) <= angles(j) + 0.33*pi
            gaindB = curve.bor(abs(possy(i) - angles(j)));
            gains(i) = gains(i) + 10^(gaindB/10);
        end
    end
end
gains = 10*log10(gains);

prism = [possy',gains'];

if plotting == 1
    polarplot(possy,gains)
end
