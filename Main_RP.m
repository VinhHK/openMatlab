%% Define Red Pitaya as TCP/IP object
IP= '169.254.182.46'; %74.85'; % Input IP of your Red Pitaya...
port = 5000;
tcpipObj=tcpip(IP, port); %tcpip not tcpclient before (This comment is due to the different matlab version that used another type code)
tcpipObj.InputBufferSize = 16384*128;
close all

%% Open connection with your Red Pitaya
fopen(tcpipObj);
tcpipObj.Terminator = 'CR/LF';
flushinput(tcpipObj);
flushoutput(tcpipObj);

fprintf(tcpipObj,'GEN:RST'); % Reset the generator parameters

%% LOAD SYSTEM REQUIRED WHEN USING THE PUREDEVICE
LoadSystem;  % load PureDevice system parameters

%% Setup Larmor frequency, LO Frequency to use,
f_Tx = HW.fLarmor; %24372832;
%pause(2) %Pause if required by multiple acquisition
f_LO = f_Tx+2000; %Set Local Oscilator frequency to be used later

Tx_Time = 66.49e-6; % Set the Tx length to be used later

%% Send channel 1 (Tx)
fprintf(tcpipObj,'SOUR1:FUNC SINE'); %Type of signal
fprintf(tcpipObj,'SOUR1:BURS:STAT BURST');      % Set type burst
Tx_f = sprintf('SOUR1:FREQ:FIX %f',f_Tx);% Set the meaning of Tx_f to be used below
fprintf(tcpipObj, Tx_f);% Set frequency of output signal
Tx_Amp = sprintf('SOUR1:VOLT %f',0.22); % Amplitude value of Source (default used: 0.22)
fprintf(tcpipObj, Tx_Amp);         % Set amplitude of output signal
Tx_length = Tx_Time * 100/5e-6; % Set the lenght of Tx
Burston = sprintf('SOUR1:BURS:NCYC %f',Tx_length); % Conversion to string for burst timing
fprintf(tcpipObj, Burston);       % Set the burst length in the Red P memory

%% Send channel 2 (LO)
fprintf(tcpipObj,'SOUR2:FUNC SIN');      % Set function of output signal {sine, square, triangle, sawu,sawd, pwm} 24381696
LO_f = sprintf('SOUR2:FREQ:FIX %f',f_LO);% Set the meaning of LO_f for the frequency of output signal
fprintf(tcpipObj, LO_f);                % Set frequency of output signal from the above LO_f
fprintf(tcpipObj,'SOUR2:VOLT 1');         % Set amplitude of output signal
fprintf(tcpipObj,'SOUR2:VOLT:OFFS 0');        %Set offset to 0
fprintf(tcpipObj,'OUTPUT2:STATE ON');        % Set output to ON

%% Receive on channel 1 (Rx)
%Set acquisition
fprintf(tcpipObj,'ACQ:RST');
fprintf(tcpipObj,'ACQ:DEC 8192');  %decimation factor
fprintf(tcpipObj,'ACQ:AVG OFF');    %disable averaging
fprintf(tcpipObj,'ACQ:TRIG:LEV 0.01');  %triger level in mV
fprintf(tcpipObj,'ACQ:SOUR1:GAIN LV'); %source and gain HV or LV

% start acquisition

fprintf(tcpipObj,'ACQ:START');
% After acquisition is started some time delay is needed in order to acquire fresh samples in to buffer
% Here we have used time delay of one second but you can calculate exact value taking in to account buffer
% length and smaling rate%

% Start measurement
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
fprintf(tcpipObj,'OUTPUT1:STATE ON'); % Set output Tx (Burst) to ON
pause(0.65); %pause (default: 0.65 with the Pure Device)
fprintf(tcpipObj,'ACQ:TRIG NOW'); % Trigger the acquisition now
% Wait for trigger PE positiv edge
% Until trigger is true wait with acquiring
% Be aware of while loop if trigger is not achieved
% Ctrl+C will stop code executing in MATLAB
% fprintf(tcpipObj,'OUTPUT1:STATE OFF'); %Set output Tx OFF if required
while 1
    trig_rsp=query(tcpipObj,'ACQ:TRIG:STAT?');
    fprintf(tcpipObj,'DIG:PIN LED3,1');
    if strcmp('TD',trig_rsp(1:2))  % Read only TD
        break
    end
    fprintf(tcpipObj,'DIG:PIN LED3,0');
end
fprintf(tcpipObj,'DIG:PIN LED3,0');

%% Read & plot
signal_str=query(tcpipObj,'ACQ:SOUR1:DATA?'); %read data from buffer
% Convert values to numbers.% First character in string is “{“
% and 2 latest are empty spaces and last is “}”.

signal_num=str2num(signal_str(1,2:length(signal_str)-3)); %str2num par default pas str2double
sample_number = length(signal_num)-1;
x_axis_raw= (0:1:sample_number);
fs=15.258e3; %sampling rate/frequency See decimation factor !!!
x_axis= x_axis_raw./fs; %sample number / Fs, for Fs see decimation factor

%% Make the FFT
FFT_signal=fft(signal_num); %Take the raw data
FFT_abs=abs(FFT_signal); %Convert the raw data to absolute
N = length(FFT_abs); % automaticaly measure the vector length
fgrid = fs*(0:(N-1))/(N);%Create the x-axis with the previous length of data
FFT_abs = FFT_abs(100:floor(N/4)); %Take only one part of the FFT
fgrid = fgrid(100:floor(N/4)); %Adapt the x-axis to the FFT taken

%% Frequency to ppm conversion
fgrid = ((fgrid*1e6)/(f_LO))-60; %Convert frequencies to PPM with the required shift to be around 0

%% Plot the signal
figure % First Figure
plot(x_axis,signal_num)
grid on
ylabel('Voltage(V)')
xlabel('Time(s)')
%hold off
figure %Second Figure
plot(fgrid,FFT_abs.^2);
ylabel('power ($V/\sqrt(Hz)$)','Interpreter','latex')
xlabel('Chemical Shift (ppm)')
grid on

%% Close connection with Red Pitaya
fclose(tcpipObj);