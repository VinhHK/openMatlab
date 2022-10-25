% Resonanter Spule?
HW.TX.Uout2PaUout(2)=HW.TX.Uout2PaUout(2)*10^(-0.2/20);
HW.TX.ChannelDef=2;         % Normaler TX HF Channal
HW.RX2TXdeadTime=5e-6;      % Totzeit des Empf�ngers vor dem Senden
HW.TX2RXdeadTime=15e-6;      % Totzeit des Empf�ngers nach dem Senden
HW.TX.BlankOffset=9e-6;     % Offset des Blankingsignals vor dem TX Puls
HW.TX.BlankPostset=1e-6;    % Zeitzusatz des Blankingsignals nach dem TX Puls
HW.RX.VGAGainDef=HW.RX.VGAGainMax/3; % Normale Empfangsverst�rkung

HW.RX.LNAGain=10^(21.7/20); % 24.8 MHz NF=1.8 dB

HW.TX.Max.PaUout(2)=min(HW.TX.Max.PaUout(2),100);

