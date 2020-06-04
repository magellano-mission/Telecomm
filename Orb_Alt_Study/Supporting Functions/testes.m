clearvars
clc
close all

sides = 4;
peak_gain = 27;
plotting = 1;

[prism] = phased_prism(sides,peak_gain,plotting);

