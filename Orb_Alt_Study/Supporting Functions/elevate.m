function [C] = elevate(a,b,c)

%% INPUTS
% a - [m] range from the transmitter to the receiver
% b - [m] distance from Mars centre to the transmitter
% c - [m] distance from Mars centre to the receiver

C = acos((a^2 + b^2 -c^2)/(2*a*b));