% LoadMySystem

HW.fLarmor = 9858000.000; HW.B0 = HW.fLarmor/(HW.Gamma.P31/2/pi);

% DC-600 external gradient amplifier
LoadGradAmp_DC600_SN_35;

% RF-100 external RF amplifier
LoadRF100_21;

% LNA - pre-amplifier and switch
LoadTxRxLNA_36;
