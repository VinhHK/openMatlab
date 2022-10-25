%% Settings for rf amplifier RF-100 with SN 21

HW.TX.ExtRFSN = 21;

HW.TX.ChannelDef = 2;         % Default TX Channel set to Tx2

% Amplifier
HW.TX.Uout2PaUout(2) = 50;       % 50x Amplification
HW.TX.Max.PaUout(2) = 101;       % maximum PaUout(2)=200 V => 201 W

HW.TX.Def.PaUout(2) = 5;        % default PaUout(2)=50 V => 100 W (10 Volt ici)

HW.TX.Max.Amplitude = [20,20]*1e-3;   % maximum B1+ in T
HW.TX.Def.Amplitude = [20,20]*1e-3;   % default B1+ in T

%% Coil 10mm
HW.RX2TXdeadTime = 1e-6;         % Receiver deadtime before TX pulse in s (old 1e-6)
HW.TX2RXdeadTime = 5e-6;        % Receiver deadtime after TX pulse in s e.g. ringing of Coil ~40 �s (Old : 5e-6)
HW.TX.BlankOffset = 800e-9;      % Unblank of RF-100 before TX pulse (old 800e-9)
HW.TX.BlankPostset = 400e-9;     % Blank of RF-100 after TX pulse (old 400e-9)
%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 1000e-9;   % Blank of receiver before TX pulse old 1000e-9
if HW.TX.DampCoil.Enable
    HW.TX.BlankPostsetAQ = 1400e-9;  % Blank of receiver after TX pulse old 1400e-9
else 
    HW.TX.BlankPostsetAQ = 2400e-9;%old 2400e-9
end
HW.TX.BlankAQ = 1;               % Switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%% RX Settings
HW.RX.VGAGainDef = HW.RX.VGAGainMax/1; % reduce VGA gain to avoid saturation

% HW.RX.LNAGain=10^(((-61.12)-(-60))/20); % 25.8757 dB gain @ 24500397.7985 MHz F=0.82651 dB (-60 dBm cal)


%%
UseRF100Switch = 1; % use switch of RF-100
LoadRF100_Cal; % new cal Uout and 6 A FET

