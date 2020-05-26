clc
clearvars

% find the pointing vectors - transmitter vector will be used to calculate
% the tranmission gain, receiver vector will be used to calculate the
% receiver gain
vec1 = [1 0.5 1];
vec2 = [-1 0 0.2];

% position vector used to calculate elevation and un-altered transmission
% and receiving angles
pos1 = [1000 500 1000];
pos2 = [5000 0 0];
rel = pos2 - pos1;

% calculate the angle between the antenna poninting directions. this may be
% useful as a check
point = atan2(norm(cross(vec1,vec2)),dot(vec1,vec2));
point = pi - point;

%calculate the final receiver incident angle
offset = atan2(norm(cross(rel,vec2)),dot(rel,vec2));
offset = pi - offset;

point = rad2deg(point);
offset = rad2deg(offset);
disp(point)
disp(offset)