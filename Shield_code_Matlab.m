%% Define Red Pitaya as TCP/IP object

close all

IP= '169.254.182.46';           % Input IP of your Red Pitaya...
port = 5000;
tcpipObj=tcpip(IP, port);
tcpipObj.InputBufferSize = 16384*16; 

%% Open connection with your Red Pitaya
fopen(tcpipObj);
tcpipObj.Terminator = 'CR/LF';
flushinput(tcpipObj);
flushoutput(tcpipObj);

% FID
% Preparations
    LoadSystem;                                     % load system parameters

    [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10, [], 1);  % find magnet frequency
     HW.fLarmor = 24373371.6266625;                              % override Larmor frequency 24.372857

    % Parameters used for timing calculations
    Seq.p90         = HW.tFlip90Def;                % duration of TX pulse
    Seq.plotSeq     = [];                           % plot sequence off

    % Sequence parameters
    Seq.tRep        = 100e-3;                       % repetition time

    % RF transmission parameters
    TX.Start        = 0;                            % start time of rf-pulse
    TX.Duration     = Seq.p90;                      % duration of rf-pulse
    TX.Frequency    = HW.fLarmor;                   % frequency of rf-pulse
    TX.Phase        = 0;                            % phase of rf-pulse

    % Acquisition parameters
    AQ.Start        = 100e-6 + Seq.p90;             % acquisition start time
    AQ.fSample      = 25e3;                         % sampling rate of AQ window
    AQ.nSamples     = 1024;                         % number of samples in AQ window
    AQ.Frequency    = HW.fLarmor;                   % frequency of AQ window
    AQ.Phase        = 0;                            % phase of AQ window
    AQ.Gain         = HW.RX(1).Amplitude2Uin / 20e-3;  % maximum input voltage
    
%% ENTETES 
f_decal=2500;
f_Tx = HW.fLarmor;
f_LO = f_Tx%+f_decal;

fprintf(tcpipObj,'GEN:RST');

%% TX
% % 
% % fprintf(tcpipObj,'DIG:PIN LED0,1'); %LED allumé pendant l'émission 
% % fprintf(tcpipObj,'SOUR1:FUNC SINE');
% % TX_send = sprintf('SOUR1:FREQ:FIX %f',f_Tx);% Set frequency of output signal
% % fprintf(tcpipObj, TX_send)
% % fprintf(tcpipObj,'SOUR1:VOLT 0.15');             % Set amplitude of output signal
% % Tx_time = HW.tFlip90Def;
% % Tx_length = Tx_time * 100/4e-6;
% % fprintf(tcpipObj,'SOUR1:BURS:STAT ON');      % Set burst mode to ON
% % Burston = sprintf('SOUR1:BURS:NCYC %f',Tx_length); % Conversion to string for burst timing
% % fprintf(tcpipObj, Burston);       % Set Tx_length to pulses of sine wave
% % 


%% LO SETTINGS
% fprintf(tcpipObj,'GEN:RST');
% 
% fprintf(tcpipObj,'SOUR1:FUNC SQUARE');       % Set function of output signal
% fprintf(tcpipObj,'SOUR2:FUNC SQUARE');       % Set function of output signal
%                                              % {sine, square, triangle, sawu,sawd, pwm}
% fprintf(tcpipObj,'SOUR1:PHAS 0');       % Phase adjust                                            
% fprintf(tcpipObj,'SOUR2:PHAS 180');       % Phase adjust
% 
% LO_1 = sprintf('SOUR1:FREQ:FIX %f',f_LO);% Set frequency of output signal
% fprintf(tcpipObj, LO_1)
% LO_2 = sprintf('SOUR2:FREQ:FIX %f',f_LO);% Set frequency of output signal
% fprintf(tcpipObj, LO_2)
% 
% fprintf(tcpipObj,'SOUR1:VOLT 1.2');         % Set amplitude of output signal
% fprintf(tcpipObj,'SOUR2:VOLT 1.2');         % Set amplitude of output signal 1.2
% fprintf(tcpipObj,'SOUR1:VOLT:OFFS 1.2');
% fprintf(tcpipObj,'SOUR2:VOLT:OFFS 1.2');

% fprintf(tcpipObj,'OUTPUT1:STATE ON');        % Set output to ON
% fprintf(tcpipObj,'OUTPUT2:STATE ON');        % Set output to ON
%% Send channel LO

fprintf(tcpipObj,'SOUR2:FUNC SINE');       % Set function of output signal {sine, square, triangle, sawu,sawd, pwm} 24381696
LO_f = sprintf('SOUR2:FREQ:FIX %f',f_LO);% Set frequency of output signal
fprintf(tcpipObj, LO_f);
fprintf(tcpipObj,'SOUR2:VOLT 1');         % Set amplitude of output signal
%fprintf(tcpipObj,'SOUR2:VOLT:OFFS 0');        
fprintf(tcpipObj,'OUTPUT2:STATE OFF');        % Set output to ON

    
%% RECEPTION

%Set acquisition
fprintf(tcpipObj,'ACQ:RST');
fprintf(tcpipObj,'ACQ:DEC 8192');  %decimation factor
fprintf(tcpipObj,'ACQ:AVG OFF')     %disable averaging
fprintf(tcpipObj,'ACQ:TRIG:LEV 0.1');  %triger level in mV
fprintf(tcpipObj,'ACQ:SOUR1:GAIN LV'); %source and gain HV or LV

%start acquisition
fprintf(tcpipObj,'ACQ:START');
%% Start measurement
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
%fprintf(tcpipObj,'OUTPUT2:STATE ON');        % Set output to ON

pause(0.7);
fprintf(tcpipObj,'ACQ:TRIG NOW'); %Trigger
%fprintf(tcpipObj,'OUTPUT2:STATE ON');        % Set output to ON

while 1
     trig_rsp=query(tcpipObj,'ACQ:TRIG:STAT?');
     fprintf(tcpipObj,'DIG:PIN LED3,1')
     if strcmp('TD',trig_rsp(1:2));  % Read only TD
     break
     end
     fprintf(tcpipObj,'DIG:PIN LED3,0')
end
     fprintf(tcpipObj,'DIG:PIN LED3,0')
     
% Read & plot
signal_str=query(tcpipObj,'ACQ:SOUR1:DATA?'); %read data from buffer
signal_num=str2num(signal_str(1,2:length(signal_str)-3));
sample_number = length(signal_num)-1;
x_axis_raw= (0:1:sample_number);
fs=15.258e3; %sampling rate/frequency See decimation factor !!!
x_axis= x_axis_raw./fs; %sample number / Fs, for Fs see decimation factor

%Faire la FFT

FFT_signal=fft(signal_num);
FFT_abs=abs(FFT_signal);
N = length(FFT_abs);
fgrid = fs*(0:(N-1))/(N);
%FFT_abs = FFT_abs(100:floor(N/2));

%fgrid = ((fgrid-f_decal)/HW.fLarmor)*(1/10e-6);

%fgrid=((((fgrid+24379745)-HW.fLarmor)/(HW.fLarmor))*10e6);
%fgrid = fgrid(100:floor(N/2));
%freq_un_ppm=(HW.fLarmor/10e6)

%ppm_coefficient_ideal = (24379745 / ((f_Tx - 24379745) / f_Tx)) / 1e6;

%fgrid = (fgrid+f_Tx) / ppm_coefficient_ideal;


%ppm_coefficient_ideal = (frequency_reception/ ((frequency_excitation - frequency_reception) / frequency_excitation)) / 1e6

%ppm_spectrum_received_ideal = frequency_spectrum_received / ppm_coefficient_ideal


%% Plot
figure
plot(x_axis,signal_num)
grid on
ylabel('Voltage(V)')
xlabel('Time(s)')
%hold off5
figure
plot(fgrid,FFT_abs);
ylabel('power / ?')
xlabel('frequencies')
grid on


% Plot results
plot_data_1D(HW, data_1D);

amplitude2rmsUcoil = data(1).Amplitude2Uin(1) / HW.RX(SeqOut.AQ(1).Device).LNAGain / sqrt(2);

hf = figure(9); clf(hf);
hax = axes(hf);
plot(hax, data(1).time_of_tRep, abs(data(1).data)*amplitude2rmsUcoil*1e6);
title(hax, 'Acquired signal');
xlabel(hax, 'time in s');
ylabel(hax, sprintf('RMS voltage at coil in %cV', char(181)));
grid(hax, 'on');

% hf = figure(10); clf(hf);
% hax = axes(hf);
% plot(hax, data.time_of_tRep, angle(data.data));
% title(hax, 'Phase of acquired signal');
% ylabel(hax, 'phase in rad');
% xlabel(hax, 'Time in s');
% grid(hax, 'on');

hf = figure(11); clf(hf);
hax = axes(hf);
% without CIC filter correction
plot(hax, ...
  data(1).f_fft1_data, ...
  abs(data(1).fft1_data)./data(1).cic_corr*amplitude2rmsUcoil*1e6);
% % with CIC filter correction
% plot(hax, data(1).f_fft1_data, abs(data(1).fft1_data).*amplitude2rmsUcoil*1e6);
xlim(hax, [data(1).f_fft1_data(1), data(1).f_fft1_data(end)])
title(hax, 'FFT of acquired signal without CIC correction');
xlabel(hax, 'frequency in Hz');
ylabel(hax, sprintf('RMS voltage at coil in %cV', char(181)));
grid(hax, 'on');

hf = figure(12); clf(hf);
hax = axes(hf);
% plot(hax, data(1).f_fft1_data-SeqOut.AQ(1).Frequency(1),abs(data(1).fft1_data)*amplitude2rmsUcoil*1e6); xlabel(hax, 'Frequency in Hz');  % offset frequency
plot(hax, (data(1).f_fft1_data/SeqOut.AQ(1).Frequency(1)-1)*1e6, ...
  abs(data(1).fft1_data)*amplitude2rmsUcoil*1e6);
title(hax, 'FFT of acquired signal with CIC correction');
xlabel(hax, 'Frequency in ppm'); % offset ppm
ylabel(hax, sprintf('RMS voltage at coil in %cV', char(181)));
grid(hax, 'on');

%% Close connection with Red Pitaya
flushinput(tcpipObj);
flushoutput(tcpipObj);
fclose(tcpipObj);