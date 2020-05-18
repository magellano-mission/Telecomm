function [out] = zfun_link_rate(frq,rang,powt,cont,G_t,G_r,dt,sit)
%% PURPOSE
%This function takes inputs relating to the link and hardware, calculates a
%link budget and then assigns a data rate based on the CNR and contware
%capabilities

%% INPUTS
%frq       - [Hz]  Carrier wave frequency
%rang      - [m]   Transmission range
%powt      - [W]   Transmitted RF power
%cont.powr - [W,W] Minimum and maximum signal power limits of receiver
%cont.M    - [-,-] Range of modulation efficiencies
%cont.P    - [-]   Polarisation of receiver and signal 'lin' or 'cir'
%cont.R    - [-]   Error encoding rate (0.5 for turbo coding)
%cont.BW   - [Hz]  Bandwidth of receiving module
%cont.Ts   - [K]   Effective noise temperature of the receiving system
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
        case 3      %Earth rain fade (NEED TO CHECK)
            mrg_atm = -3;
    end
end
if frq > 4e9 && frq < 12e9 %X band frequency range
    switch sit
        case 1      %Martian atmospheric losses
            mrg_atm = - 1.15;
        case 2      %No losses (S/C to S/C)
            mrg_atm = 0;
        case 3      %Earth X (NEED TO CHECK)
            mrg_atm = -1;  %not including rain fade
    end
end
if frq > 30e9 && frq < 40e9 %Ka band frequency range
    switch sit
        case 1      %Martian atmospheric losses
            mrg_atm = - 3.35;
        case 2      %No losses (S/C to S/C)
            mrg_atm = 0;
        case 3      %Earth Ka (NEED TO CHECK)
            mrg_atm = -2;   %not including rain fade
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
if PdB_r >= cont.powr(1) && PdB_r <= cont.powr(2)
    %Calculating CNR
    n0   = physconst('Boltzmann')*(cont.Ts);    %[W/Hz] Noise power spectral density
    n0_W = n0*cont.BW;                          %[W] Noise power
    n0_dB = 10*log10(n0_W);                     %[dBW] Noise Power
    CNRdB = PdB_r - n0_dB;                      %[dB] Carrier to noise power ratio
    CNR = 10^(CNRdB/10);                        %[-] converting CNR to linear ratio
    fb = cont.symmax * log2(cont.M);            %[bits/s] maximum bit rate possible
    EbNo = CNR * cont.BW / fb;                  %[-] Bit energy to noise power spectral density ratio
    ENdB = 10*log10(EbNo);                      %[-] Converted to log form

    % Finding the acceptable Eb/No ratio for the particular coding type.
    % Eb/No [dB] limits in order to keep BER < 10^-6, from Satellite
    % Communications textbook
    % MPSK assumed based on hardware options examined
    switch cont.M
            case 2  
                ENdB_lim = 10.6;
            case 4 
                ENdB_lim = 10.6;
            case 8
                ENdB_lim = 14.0;
            case 16
                ENdB_lim = 18.3;    
    end
    
    % Check if the Eb/No is above the acceptable limit
    if ENdB >= ENdB_lim
        rate = fb;
    else            %Introduce error coding to reduce Eb/No limit
       ENdB_lim = ENdB_lim - cont.cgain;
       if ENdB >= ENdB_lim
           rate = fb * cont.R;
       else         %Iteratively reduce the bit rate
           rate = fb * cont.R;
           while ENdB <= ENdB_lim 
               rate = rate/2;
               ENdB = ENdB + 3;
           end      %Check if final bit rate is compatible with the hardware
           if rate <= cont.symmin * log2(cont.M)
               rate = 0;
           end
       end
    end
else
    %Otherwise transmission is not possible
    CNRdB = NaN;
    rate = 0;
end
 
%% OUTPUT
%Saving variables to an output structure
out.PdB   = PdB_r;                %[dBW] Received power
out.P_W   = PW_r;                 %[W]   Received power raw
out.dE    = out.P_W * dt;         %[J]   Received energy in time step
out.CNRdB = CNRdB;                %[dB]  log10 of CNR power ratio
out.datrat = rate;                %[bits/s] rate of useful information transfer
out.data = rate*dt;               %[bit] volume of data transferred in time step
