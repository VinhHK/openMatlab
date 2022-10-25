HW.Grad.ExtGradSN = 31;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[0.0001372688,-0.0001198268,0.0001923585,0.0007341912]; % Offset Current 15-Mar-2021 14:02:42
 
HW.Grad.PaUin2PaIout=([0.3329833,0.3319632,0.3323263,0.3326932]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 15-Mar-2021 14:02:42

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time in s
HW.Grad.tEC=50e-6;                                      % Setting time;

switch HW.UserName
  case 'm175_probe_H1_10'
    % Magnet 175 - 10 mm RF coil - 10 mm Gradsystem V2    
    HW.Grad.SystemTimeDelay(1:3) = [2.29136e-05   3.1116e-05  2.66376e-05]; % Time delay of grad amp

  case 'm175_probe_H1_5'
    % Magnet 175 - 5 mm RF coil - 10 mm Gradsystem V2
    HW.Grad.SystemTimeDelay(1:3) = [2.36592e-05  3.24192e-05  3.49088e-05]; % Time delay of grad amp

  otherwise
    % Magnet 175 - 10 mm RF coil - 10 mm Gradsystem V2
    HW.Grad.SystemTimeDelay(1:3) = [2.29136e-05   3.1116e-05  2.66376e-05]; % Time delay of grad amp

end

HW.Grad.MaxAmpSlice=0.1;                                % Max Grad Amp ?

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600
