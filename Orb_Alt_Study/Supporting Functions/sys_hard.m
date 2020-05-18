function [hard] = sys_hard(antt,antr,custt,custr,cont,T_amb,point)
%% PURPOSE
%This function collects the relevant options and information about the
% possible hardware selection options

%% ANTENNA OPTIONS
%Antenna gain curves vs boresight angle and angle limit information
[PAX] = MES_PAX_8490(0);
[UHFB] = MRO_UHF_401_6(0);   %[dBi] MRO / Curiosity UHF antenna curve
[UHFI] = Improved_UHF_401_6(0);   %[dBi] More directional UHF antenna
[UHF2I] = Improved2_UHF_401_6(0);   %[dBi] More directional UHF antenna
[LGAX] = MRO_LGA_7183(0);    %[dBi] MRL X band LGA  curve
[MGAX] = MSL_MGA_7183(0);    %[dBi] MSL X band MGA curve
[RHGAX] = Cur_HGA_7183(0);    %[dBi] Curiosity X band HGA curve
[OHGAX] = MRO_HGA_7183(0);  %[dBi] MRO X band HGA curve
[DSN34X] = DSN_34m_7183(0); %[dBi] DSN 34m dish gain (curve not available)
[OHGAKa] = MRO_Ka_34k(0);
[DSN34Ka] = DSN_34m_34k(0);

%% ANTENNA SELECTION
if isstruct(custt) == 0
    switch antt
        case 'PAX'
            hard.gaint = PAX;
        case 'UHFB'
            hard.gaint = UHFB;
        case 'UHFI'
            hard.gaint = UHFI;
        case 'UHF2I'
            hard.gaint = UHF2I;
        case 'LGAX'
            hard.gaint = LGAX;
        case 'MGAX'
            hard.gaint = MGAX;
        case 'RHGAX'
            hard.gaint = RHGAX;
        case 'OHGAX'
            hard.gaint = OHGAX;
        case 'DSN34X'
            hard.gaint = DSN34X;
        case 'OHGAKa'
            hard.gaint = OHGAKa;
        case 'DSN34Ka'
            hard.gaint = DSN34Ka;
    end
else
    [hard.gaint] = basic_antenna(custt.type,custt.gain_peak, custt.HPBW, custt.plotting);
end
if isstruct(custr) == 0
    switch antr
        case 'PAX'
            hard.gainr = PAX;
        case 'UHFB'
            hard.gainr = UHFB;
        case 'UHFI'
            hard.gainr = UHFI;
        case 'UHF2I'
            hard.gainr = UHF2I;
        case 'LGAX'
            hard.gainr = LGAX;
        case 'MGAX'
            hard.gainr = MGAX;
        case 'RHGAX'
            hard.gainr = RHGAX;
        case 'OHGAX'
            hard.gainr = OHGAX;
        case 'DSN34X'
            hard.gainr = DSN34X;
        case 'OHGAKa'
            hard.gainr = OHGAKa;
        case 'DSN34Ka'
            hard.gainr = DSN34Ka;
    end
else
    [hard.gainr] = basic_antenna(custr.type,custr.gain_peak, custr.HPBW, custr.plotting);
end

%% ANTENNA POINTING
hard.point = point;
%% CONTROLLER OPTIONS
%Electra specifications for UHF baseline
elec.powr = [-170 -100];    %[dBW] Acceptable signal power range to feed Electra
elec.symmax = 2.048e6;      %[symbols/sec]
elec.symmin = 1000;         %[symbols/sec]
elec.M = 2;                 %[bit/sym] Max modulation efficiency
elec.R = 0.5;               %[-] Error coding efficiency
elec.cgain = 6.6;            %[dB] Error coding gain
elec.BW = elec.symmax/1.7;     %[Hz] Electra bandwidth with ballpark correction factor
elec.P = 'circ';            %[-]Polarisation of signal, 'circ' or 'lin'
elec.F = 3.9;                         %[dB] Noise Factor for Electra (half-duplex)
elec.Ts = T_amb*(10^(elec.F/10)-1);   %[K] System effective noise temperature

%SDST specifications for X and Ka band baseline
sdst.powr = [-188 -100];     %[dBW] Acceptable signal power range to feed SDST
sdst.symmax = 15e6;          %[symbols/sec]
sdst.symmin = 1000;          %[symbols/sec]
sdst.M = 4;                  %[bit/sym] Max modulation efficiency (BPS,QPS,etc)
sdst.R = 0.5;                %[-] Error coding efficiency
sdst.cgain = 6.6;            %[dB] Error coding gain
sdst.BW = sdst.symmax / 1.7;    %[Hz] SDST bandwidth with ballpark correction factor
sdst.P = 'circ';             %[-]Polarisation of signal, 'circ' or 'lin'
sdst.F = 2.5;                         %[dB] Noise Factor for SDST
sdst.Ts = T_amb*(10^(sdst.F/10)-1);   %[K] System effective noise temperature

%DSN specifications for X and Ka band baseline
dsn.powr = [-300 -50];
dsn.symmax = 15e6;
dsn.symmin = 1000;
dsn.M = 4;
dsn.R = 0.5;
dsn.cgain = 10; %6.6;            
dsn.BW = dsn.symmax / 1.7;
dsn.P = 'circ';
dsn.Ts = 20;

%% CONTROLLER SELECTION
switch cont
    case 'elec'
        hard.cont = elec;
    case 'sdst'
        hard.cont = sdst;
    case 'dsn'
        hard.cont = dsn;
end