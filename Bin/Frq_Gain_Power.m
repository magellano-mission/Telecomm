clearvars
clc
tic

D_tra = 0.1;     %[dB] transmitter antenna diameter
D_rec = 5;     %[dB] receiver antenna diameter
dist  = linspace(200,10e3,50);  %[km] transmission distance

%freqGHz = [0.3 1 2 4 8 12 27 40];    %[GHz] range of frequencies
%names = ["UHF Lower" "UHF Upper" "S Lower" "S Upper" "X Lower" ...
    %"X Upper" "Ka Lower" "Ka Upper"];
freqGHz = linspace(0.3,40,50);

%Unit Conversions
dist_m = 1000*dist;                             %[m] converting km to m
freqHz = 1e9 * freqGHz;                         %[Hz] converting GHz to Hz
lambda = physconst('LightSpeed')./freqHz;       %[m] calculating wavelength
G_tra = 20*log10(pi*D_tra./lambda);             %[dB] transmission gain
G_rec = 20*log10(pi*D_rec./lambda);             %[dB] receiver gain
G_sys = G_tra + G_rec;                          %[db} combined antenna gain

%Calculations
ratio = struct;
i = 1;
hold on
for j = 1:length(dist_m)
    for k = 1:length(lambda)
    path_loss = 20*log10(lambda(k)/(4*pi*dist(j)));
    ratio.s(j,k) = G_sys(k) + path_loss;
    end
end
surf(freqGHz,dist,ratio.s)
xlabel('Frequency [GHz]')
ylabel('Distance [km]')
zlabel('Power Ratio')
toc
hold off