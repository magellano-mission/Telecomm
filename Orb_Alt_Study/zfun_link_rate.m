function [out] = zfun_link_rate(frq,dis,powt,powr,G_t,G_r,T_s,BW,dt,code)
%
%polarisation loss
%atmosphere losses

%% INPUTS and OUTPUTS
%frq -  [Hz]  Carrier wave frequency
%dis -  [m]   Transmission distance
%powt - [W]   Transmitted RF power
%powr - [W]   Minimum signal power at receiver
%G_t -  [dBi] Transmitting antenna gain
%G_r    [dBi] Receiving antenna gain
%T_s    [K]   Effective noise temperature of the receiving system
%pol    [-]   Polarisation of receiver and signal 'lin' or 'cir'
%R      [-]   Error encoding rate (0.5 for turbo coding)
%dt     [s]   Time step

%%
lambda = physconst('LightSpeed')/frq;   %[m] wavelength of signal
k = physconst('Boltzmann');

PdB_t = 10*log10(powt);      %[dBW] Transmitted signal power
pathdB_L = -20*log10(4*pi*dis/lambda);  %[dB] Propagation Loss
PdB_r = PdB_t + G_t + pathdB_L + G_r;   %[dBW] Received signal power
PW_r =  powt*10^(PdB_r/10);             %[W] Received signal power

%Checking to see if received power greater than hardware minimum
if PdB_r >= powr(1) && PdB_r <= powr(2)
    %Calculating CNR for the given bandwidth
    CNR = PW_r/(k*T_s*BW);
    CNRdB = 10*log10(CNR);
else
    %Transmission not possible
    CNR = NaN;
end

%Finding the best possible symbol and bit rates for the given CNR
if CNRdB < 5
    M = 0; 
    R =0;
end
if CNRdB >= 5 && CNRdB < 12
    M = code.M(1);
    R = code.R(1);
end
if CNRdB >= 12 && CNRdB < 15
    M = code.M(1);
    R = code.R(2);
end
if CNRdB >= 15
    M = code.M(2);
    R = code.R(2);
end

rate = R * M * BW;
data = rate*dt;
 
%Saving variables to an output structure for use in calling script
out.PdB   = PdB_r;                %[dBW] Received power 
out.P_W   = PW_r;                 %[W]   Received power
out.E     = out.P_W * dt;         %[J]   Received energy in time step
out.CNR   = CNR;                  %[-]   Raw CNR power ratio
out.CNRdB = CNRdB;        %[dB]  log of CNR power ratio
out.rate = rate;          %[bit/s] of useful data
out.data = data;          %[bit] of useful data in the time period

