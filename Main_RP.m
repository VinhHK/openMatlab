%% Define Red Pitaya as TCP/IP object
IP= '169.254.182.46'; %74.85'; % Input IP of your Red Pitaya...
port = 5000;
RP_Client=tcpclient(IP, port); %tcpip not tcpclient before (This comment is due to the different matlab version that used another type code)
RP_Client.ByteOrder = "big-endian";
close all

RP_Client = tcpclient(IP,port);
configureTerminator(RP_Client,"CR/LF");


%% Open connection with your Red Pitaya
fopen(RP_Client);

flush(RP_Client,"input")
flush(RP_Client,"output")

writeline(RP_Client,'GEN:RST'); % Reset the generator parameters

%% LOAD SYSTEM REQUIRED WHEN USING THE PUREDEVICE
LoadSystem;  % load PureDevice system parameters

%% Setup Larmor frequency, LO Frequency to use,
f_Tx = HW.fLarmor; %24372832;
%pause(2) %Pause if required by multiple acquisition
f_LO = f_Tx+2000; %Set Local Oscilator frequency to be used later

Tx_Time = 66.49e-6; % Set the Tx length to be used later

%% Send channel 1 (Tx)
writeline(RP_Client,'SOUR1:FUNC SINE'); %Type of signal
writeline(RP_Client,'SOUR1:BURS:STAT BURST');      % Set type burst
Tx_f = sprintf('SOUR1:FREQ:FIX %f',f_Tx);% Set the meaning of Tx_f to be used below
writeline(RP_Client, Tx_f);% Set frequency of output signal
Tx_Amp = sprintf('SOUR1:VOLT %f',0.22); % Amplitude value of Source (default used: 0.22)
writeline(RP_Client, Tx_Amp);         % Set amplitude of output signal
Tx_length = Tx_Time * 100/5e-6; % Set the lenght of Tx
Burston = sprintf('SOUR1:BURS:NCYC %f',Tx_length); % Conversion to string for burst timing
writeline(RP_Client, Burston);       % Set the burst length in the Red P memory

%% Send channel 2 (LO)
writeline(RP_Client,'SOUR2:FUNC SIN');      % Set function of output signal {sine, square, triangle, sawu,sawd, pwm} 24381696
LO_f = sprintf('SOUR2:FREQ:FIX %f',f_LO);% Set the meaning of LO_f for the frequency of output signal
writeline(RP_Client, LO_f);                % Set frequency of output signal from the above LO_f
writeline(RP_Client,'SOUR2:VOLT 1');         % Set amplitude of output signal
writeline(RP_Client,'SOUR2:VOLT:OFFS 0');        %Set offset to 0
writeline(RP_Client,'OUTPUT2:STATE ON');        % Set output to ON

%% Receive on channel 1 (Rx)
%Set acquisition
writeline(RP_Client,'ACQ:RST');
writeline(RP_Client,'ACQ:DEC 8192');  %decimation factor
writeline(RP_Client,'ACQ:AVG OFF');    %disable averaging
writeline(RP_Client,'ACQ:TRIG:LEV 0.01');  %triger level in mV
writeline(RP_Client,'ACQ:SOUR1:GAIN LV'); %source and gain HV or LV

% start acquisition

writeline(RP_Client,'ACQ:START');
% After acquisition is started some time delay is needed in order to acquire fresh samples in to buffer
% Here we have used time delay of one second but you can calculate exact value taking in to account buffer
% length and smaling rate%

% Start measurement
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
writeline(RP_Client,'OUTPUT1:STATE ON'); % Set output Tx (Burst) to ON
pause(0.65); %pause (default: 0.65 with the Pure Device)
writeline(RP_Client,'ACQ:TRIG NOW'); % Trigger the acquisition now
% Wait for trigger PE positiv edge
% Until trigger is true wait with acquiring
% Be aware of while loop if trigger is not achieved
% Ctrl+C will stop code executing in MATLAB
% writeline(RP_Client,'OUTPUT1:STATE OFF'); %Set output Tx OFF if required
while 1
    trig_rsp=query(RP_Client,'ACQ:TRIG:STAT?');
    writeline(RP_Client,'DIG:PIN LED3,1');
    if strcmp('TD',trig_rsp(1:2))  % Read only TD
        break
    end
    writeline(RP_Client,'DIG:PIN LED3,0');
end
writeline(RP_Client,'DIG:PIN LED3,0');

%% Read & plot
signal_str=query(RP_Client,'ACQ:SOUR1:DATA?'); %read data from buffer
% Convert values to numbers.% First character in string is “{“
% and 2 latest are empty spaces and last is “}”.

signal_num=str2num(signal_str(1,2:length(signal_str)-3)); %str2num par default pas str2double
sample_number = length(signal_num)-1;
x_axis_raw= (0:1:sample_number);
fs=15.258e3; %sampling rate/frequency See decimation factor !!!
x_axis= x_axis_raw./fs; %sample number / Fs, for Fs see decimation factor

%% Create the imaginary part 

Fid_IQ = hilbert(signal_num); % Using hilbert transform we have the 90 Phase signal (Quadrature) en the real part


%% Make the FFT
FFT_signal_IQ=fft(Fid_IQ); %Take the raw data
FFT_abs=abs(FFT_signal_IQ); %Convert the raw data to absolute
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
clear(RP_Client);