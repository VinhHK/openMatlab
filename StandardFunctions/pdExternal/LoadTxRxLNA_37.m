% Load settings for LNA #37 24.15 MHz
HW.RX.LnaSN = 37;

HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.3/20); % -0.3 dB Gain (TX to Coil) @ 24.15 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 10e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 9e-6; % >7�s         % blank of transmit before TX pulse
HW.TX.BlankPostset = 5e-6;            % blank of transmit after TX pulse

HW.TX.BlankOffsetAQ = 2e-6;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 5e-6;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain
HW.RX.LNAGain=10^(((-44.0412)-(-67.6))/20); % 23.5588 dB gain @ 24150350.4346 MHz F=0.72833 dB (-67.6 dBm cal)
 
HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 10);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 4);  % def transmit voltage