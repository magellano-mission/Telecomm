function [out] = zfun_link_rate(frq,rang,powt,powr,G_t,G_r,T_s,BW,dt,code)
%% TO ADD LIST
%polarisation loss
%atmosphere losses

%% INPUTS
%frq -  [Hz]  Carrier wave frequency
%rang - [m]   Transmission range
%powt - [W]   Transmitted RF power
%powr - [W,W] Minimum and maximum signal power limits of receiver
%G_t -  [dBi] Transmitting antenna gain (elevation included)
%G_r    [dBi] Receiving antenna gain (elevation included)
%T_s    [K]   Effective noise temperature of the receiving system
%code.M [-,-] Range of modulation efficiencies
%code.P [-]   Polarisation of receiver and signal 'lin' or 'cir'
%code.R [-]   Error encoding rate (0.5 for turbo coding)
%dt     [s]   Time step

%% CALCULATIONS
lambda = physconst('LightSpeed')/frq;    %[m] wavelength of signal

%Link budget for received power
PdB_t = 10*log10(powt);                  %[dBW] Transmitted signal power
pathdB_L = -20*log10(4*pi*rang/lambda);  %[dB] Propagation Loss
PdB_r = PdB_t + G_t + pathdB_L + G_r;    %[dBW] Received signal power
PW_r =  powt*10^(PdB_r/10);              %[W] Received signal power

%Checking if received signal power is within hardware limits
if PdB_r >= powr(1) && PdB_r <= powr(2)
    %Calculating CNR for the given bandwidth
    CNR = PW_r/(physconst('Boltzmann')*T_s*BW);
    CNRdB = 10*log10(CNR);
else
    %Transmission not possible
    CNR = NaN;
    CNRdB = NaN;
end

%Finding the best possible symbol and bit rates for the given CNR in dB
% CNRdB values from table in Satellite Communication book

% CURRENTLY LIMITED TO QPSK, REVISIT THIS SECTION IF HIGHER RATES 
if CNRdB < 5 || isnan(CNRdB)
    M = 0;          %no data transmission
    R =0;           %no data transmission
end
if CNRdB >= 5 && CNRdB < 12
    M = 2;          %BPSK
    R = code.R;     %turbo coding applied for error correction
end
if CNRdB >= 12 && CNRdB < 15
    M = 2;          %BPSK
    R = 1;          %no error correction coding
end
if CNRdB >= 15
    if code.M >= 4
    M = 4;          %QPSK
    R = 1;          %no error correction coding
    else
        M= code.M;
        R = 1;
    end
end
if CNRdB >= 19
    if code.M >= 8
    M = 8;          %8PSK
    R = 1;          %no error correction coding
    else
        M= code.M;
        R = 1;
    end
end

rate = R * M * BW;  %[bits/s] rate of useful information transfer
data = rate*dt;     %[bits] volume of data transferred in time step
 
%% OUTPUT
%Saving variables to an output structure
out.PdB   = PdB_r;                %[dBW] Received power
out.P_W   = PW_r;                 %[W]   Received power raw
out.dE    = out.P_W * dt;         %[J]   Received energy in time step
out.CNR   = CNR;                  %[-]   Raw CNR power ratio
out.CNRdB = CNRdB;                %[dB]  log10 of CNR power ratio
out.datrat = rate;                  %[bits/s] rate of useful information transfer
out.data = data;                  %[bit] volume of data transferred in time step
