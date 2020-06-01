function [C] = elevation(a,b,c)

C = acos((a^2 + b^2 -c^2)/(2*a*b));