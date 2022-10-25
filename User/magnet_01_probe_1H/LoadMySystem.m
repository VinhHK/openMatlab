% LoadMySystem

HW.fLarmor = 24371000.000; %HW.B0 = HW.fLarmor/(HW.Gamma.H1/2/pi);
%HW.TX(1).PaUout2Amplitude = [0.001, 10.462705]*1e-6;  % 2022-05-10T14:36:43 (tFlip90 = 2.373 µs @ 10.000 V) from 1d Spin 0.00015, 247.433717
% DC-600 external gradient amplifier
LoadGradAmp_DC600_SN_35;    

% RF-100 external RF amplifier
LoadRF100_21;

% LNA - pre-amplifier and switch
%LoadTxRxLNA_36;
