% LoadMySystem

HW.fLarmor = 24380000.000; HW.B0 = HW.fLarmor/(HW.Gamma.H1/2/pi);

% RF-100 external RF amplifier
LoadRF100_13;

% LNA - pre-amplifier and switch
LoadTxRxLNA_24;
