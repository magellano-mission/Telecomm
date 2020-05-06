function [out] = zfun_link_rate(frq,rang,powt,hard,G_t,G_r,BW,dt,sit)
%% PURPOSE
%This function takes inputs relating to the link and hardware, calculates a
%link budget and then assigns a data rate based on the CNR and hardware
%capabilities

%% INPUTS
%frq       - [Hz]  Carrier wave frequency
%rang      - [m]   Transmission range
%powt      - [W]   Transmitted RF power
%hard.powr - [W,W] Minimum and maximum signal power limits of receiver
%hard.M    - [-,-] Range of modulation efficiencies
%hard.P    - [-]   Polarisation of receiver and signal 'lin' or 'cir'
%hard.R    - [-]   Error encoding rate (0.5 for turbo coding)
%hard.Ts   - [K]   Effective noise temperature of the receiving system
%G_t       - [dBi] Transmitting antenna gain (elevation included)
%G_r       - [dBi] Receiving antenna gain (elevation included)
%BW        - [Hz]  Received bandwidth being considered 
%dt        - [s]   Time step
%sit       - [-]   1 - Mars atmosphere, 2 - No atmosphere, 3 - Earth atmosphere

%% ATMOSPHERIC LOSSES
if frq < 500e6      %UHF frequency range
    switch sit
        case 1      %Martian atmospheric losses
            mrg_atm = -0.5;
        case 2      %No losses (S/C to S/C)
            mrg_atm = 0;
        case 3      %Earth atmospheric losses (NEED TO CHECK)
            mrg_atm = -0.5;
    end
end
if frq > 4e9 && frq < 12e9 %X band frequency range
    switch sit
        case 1      %Martian atmospheric losses
            mrg_atm = - 1.0;
        case 2      %No losses (S/C to S/C)
            mrg_atm = 0;
        case 3      %Earth atmospheric losses (NEED TO CHECK)
            mrg_atm = -3;
    end
end

%% CALCULATIONS
lambda = physconst('LightSpeed')/frq;    %[m] wavelength of signal

%Link budget margins
mrg_op = -1;          %[dB] typical operating margin applied
mrg_bo = -1;          %[dB] typical transponder back-off
mrg_pol = 0;          %[dB] polarisation loss

PdB_t = 10*log10(powt);                  %[dBW] Transmitted signal power
pathdB_L = -20*log10(4*pi*rang/lambda);  %[dB] Propagation Loss

%Link budget for received power
PdB_r = PdB_t + G_t + pathdB_L + G_r + mrg_op + mrg_bo + mrg_atm + mrg_pol;
PW_r =  powt*10^(PdB_r/10);              %[W] Received signal power

%Checking if received signal power is within hardware limits
if PdB_r >= hard.powr(1) && PdB_r <= hard.powr(2)
    %Calculating CNR for the given bandwidth
    n0_W = (physconst('Boltzmann')*(hard.Ts)*BW);
    n0_dB = 10*log10(n0_W);
    CNRdB = PdB_r - n0_dB;
else
    %Otherwise transmission is not possible
    CNRdB = NaN;
end

%Finding the best possible symbol and bit rates for the given CNR in dB
% CNRdB values from table in Satellite Communication book

if CNRdB < 3 || isnan(CNRdB)
    M = 0;          %no data transmission
    R =0;           %no data transmission
end
if CNRdB >= 3 && CNRdB < 11
    M = 2;          %BPSK
    R = hard.R;     %turbo coding applied for error correction
end
if CNRdB >= 11 && CNRdB < 14
    M = 2;          %BPSK
    R = 1;          %no error correction coding
end
if CNRdB >= 14 && CNRdB < 19
    if hard.M >= 4
    M = 4;          %QPSK
    R = 1;          %no error correction coding
    else
        M= hard.M;
        R = 1;
    end
end
if CNRdB >= 19 && CNRdB < 24
    if hard.M >= 8
    M = 8;          %8PSK
    R = 1;          %no error correction coding
    else
        M= hard.M;
        R = 1;
    end
end
if CNRdB >= 24
    if hard.M >= 16
    M = 16;         %16PSK (sdst can do 100 Mb/s)
    R = 1;          %no error correction coding
    else
        M= hard.M;
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
out.CNRdB = CNRdB;                %[dB]  log10 of CNR power ratio
out.datrat = rate;                %[bits/s] rate of useful information transfer
out.data = data;                  %[bit] volume of data transferred in time step
