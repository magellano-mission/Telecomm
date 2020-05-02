function [output] = zfun_link_budget(frq,BW,dist,pow_t,dia_t,e_a,atm_t, ...
    dia_r,k_r)
%This function will take input variables about the link in question and
%populate the nominal link budget.
%This code should be compared with a couple of online calculators and
%reference papers such as MRO to ensure it is correct.
%BER - bit error rate
%SNR - Signal to noise ratio
%CNR - carrier to noise ratio

lambda = physconst('LightSpeed')/frq;   %[m] wavelength of signal

PdB_t = 10*log10(pow_t);      %[dBW] Transmitted signal power
G_t   = e_a*(pi*dia_t/lambda)^2;     %[-] Transmitting antenna gain
GdB_t = 10*log10(G_t);        %[dBi] Gain relative to isotropic transmitter

pathdB_L = -20*log10(4*pi*dist/lambda);  %[dB] Propagation Loss

G_r = k_r*(pi*dia_r/lambda)^2;       %[-] Receiving antenna gain
GdB_r = 10*log10(G_r);        %[dBi] Gain relative to isotropic receiver
PdB_r = PdB_t + GdB_t + pathdB_L + GdB_r;   %[dBW] Received signal power
%Noise temperature etc.  %physconst('Boltzmann'); 

%Saving variables to an output structure for use in calling script
output.PdB_t = PdB_t;
output.GdB_t = GdB_t;
output.pathdB_L = pathdB_L;
output.GdB_r = GdB_r;
output.PdB_r = PdB_r;


