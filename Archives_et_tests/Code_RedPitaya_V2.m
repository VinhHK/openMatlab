%% Define Red Pitaya as TCP/IP object

IP= '169.254.182.46';           % Input IP of your Red Pitaya...
port = 5000;
tcpipObj=tcpip(IP, port); %tcpip pas tcpclient avant
tcpipObj.InputBufferSize = 16384*128; 
close all
%% Open connection with your Red Pitaya
fopen(tcpipObj);
tcpipObj.Terminator = 'CR/LF';
flushinput(tcpipObj);
flushoutput(tcpipObj);
%% Send a sine signal burst Tx 
% 

% for k=1:100:50000;
%     {

fprintf(tcpipObj,'GEN:RST');
% fprintf(tcpipObj,'DIG:PIN LED0,1'); %LED allumé pendant l'émission 
% fprintf(tcpipObj,'SOUR2:FUNC SINE');
% % frequ=21614740;
% % frequ_mix=21100000;
% 
% fprintf(tcpipObj,'SOUR2:VOLT 1');             % Set amplitude of output signal
% 
% Tx_time = 40.652e-6;
% 
% Tx_length = Tx_time * 100/4e-6;
% 
% fprintf(tcpipObj,'SOUR2:BURS:STAT ON');      % Set burst mode to ON
% Burston = sprintf('SOUR2:BURS:NCYC %f',Tx_length); % Conversion to string for burst timing
% fprintf(tcpipObj, Burston);       % Set Tx_length to pulses of sine wave
% %fprintf(tcpipObj,'SOUR1:BURS:NOR 1000');    % Infinity number of sine wave pulses
% %fprintf(tcpipObj,'SOUR1:BURS:INT:PER 10');  % Set time of burst period in microseconds = 5 * 1/Frequency * 1000000


%  fprintf(tcpipObj, 'DIG:PIN LED0,0'); %LED eteinte en fin d'émission
 
% %% FID
% % Preparations
    LoadSystem;                                     % load system parameters
% 
%     [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10, [], 1);  % find magnet frequency
%     HW.fLarmor = 10e6;                              % override Larmor frequency
% 
%  %   Parameters used for timing calculations
%     Seq.p90         = HW.tFlip90Def;                % duration of TX pulse
%     Seq.plotSeq     = [];                           % plot sequence off
% 
%  %   Sequence parameters
%     Seq.tRep        = 100e-3;                       % repetition time
% 
%  %   RF transmission parameters
%     TX.Start        = 0;                            % start time of rf-pulse
%     TX.Duration     = 1e-6 ;%Seq.p90;                      % duration of rf-pulse
%     TX.Frequency    = HW.fLarmor;                   % frequency of rf-pulse
%     TX.Phase        = 0;                            % phase of rf-pulse
% 
% %    Acquisition parameters
%     AQ.Start        = 100e-6 + Seq.p90;             % acquisition start time 100e-6
%     AQ.fSample      = 25e3;                         % sampling rate of AQ window
%     AQ.nSamples     = 1024;                         % number of samples in AQ window
%     AQ.Frequency    = HW.fLarmor;                   % frequency of AQ window
%     AQ.Phase        = 0;                            % phase of AQ window
%     AQ.Gain         = HW.RX(1).Amplitude2Uin / 20e-3;  % maximum input voltage
%     

f_Tx = 24371296;%24375780
f_LO = f_Tx;

% %% Send channel 1
            fprintf(tcpipObj,'S5OUR1:FUNC SINE');
            Tx_f = sprintf('SOUR1:FREQ:FIX %f',f_Tx);% Set frequency of output signal
            fprintf(tcpipObj, Tx_f);
            Tx_Amp = sprintf('SOUR1:VOLT %f',0.22); % Amplitude value of Source 2 42
            fprintf(tcpipObj, Tx_Amp);         % Set amplitude of output signal
            %writeline(app.rp,'SOUR1:VOLT:OFFS 0');
            Tx_length = 80e-6 * 100/4.05e-6;
            fprintf(tcpipObj,'SOUR1:BURS:STAT ON');      % Set burst mode to ON
            Burston = sprintf('SOUR1:BURS:NCYC %f',Tx_length); % Conversion to string for burst timing
            fprintf(tcpipObj, Burston);       % Set Tx_length to pulses of sine wave
%% Send channel LO

fprintf(tcpipObj,'SOUR2:FUNC SINE');      % Set function of output signal {sine, square, triangle, sawu,sawd, pwm} 24381696
LO_f = sprintf('SOUR2:FREQ:FIX %f',f_LO);% Set frequency of output signal
fprintf(tcpipObj, LO_f);
fprintf(tcpipObj,'SOUR2:VOLT 1');         % Set amplitude of output signal
fprintf(tcpipObj,'SOUR2:VOLT:OFFS 0');        
fprintf(tcpipObj,'OUTPUT2:STATE OFF');        % Set output to ON

%% Receive on channel 1
%Set acquisition
fprintf(tcpipObj,'ACQ:RST');
fprintf(tcpipObj,'ACQ:DEC 8192');  %decimation factor
fprintf(tcpipObj,'ACQ:AVG OFF');    %disable averaging
fprintf(tcpipObj,'ACQ:TRIG:LEV 0.1');  %triger level in mV
fprintf(tcpipObj,'ACQ:SOUR1:GAIN LV'); %source and gain HV or LV
%fprintf(tcpipObj,'ACQ:DATA:FORMAT BIN');
%fprintf(tcpipObj,'ACQ:DATA:UNITS VOLTS');
% Set trigger delay to 0 samples
% 0 samples delay set trigger to center of the buffer
% Signal on your graph will have trigger in the center (symmetrical)
% Samples from left to the center are samples before trigger
% Samples from center to the right are samples after trigger
%fprintf(tcpipObj,'ACQ:TRIG:DLY:NS 40000000000'); %delay en ns NS 17000 or 8192
%start acquisition

%fprintf(tcpipObj,'SOUR1:TRIG:IMM');          % Set generator trigger to immediately

fprintf(tcpipObj,'ACQ:START');
% After acquisition is started some time delay is needed in order to acquire fresh samples in to buffer
% Here we have used time delay of one second but you can calculate exact value taking in to account buffer
% length and smaling rate%

% Start measurement
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
fprintf(tcpipObj,'OUTPUT1:STATE ON');        % Set output to ON
pause(0.5);
fprintf(tcpipObj,'ACQ:TRIG NOW'); %Tx
% fprintf(tcpipObj,'OUTPUT1:STATE OFF');
% Wait for trigger PE positiv edge
% Until trigger is true wait with acquiring
% Be aware of while loop if trigger is not achieved
% Ctrl+C will stop code executing in MATLAB
% fprintf(tcpipObj,'OUTPUT1:STATE OFF');
while 1
     trig_rsp=query(tcpipObj,'ACQ:TRIG:STAT?');
     fprintf(tcpipObj,'DIG:PIN LED3,1');
     if strcmp('TD',trig_rsp(1:2))  % Read only TD
     break
     end
     fprintf(tcpipObj,'DIG:PIN LED3,0');
end
     fprintf(tcpipObj,'DIG:PIN LED3,0');
     
% Read & plot

signal_str=query(tcpipObj,'ACQ:SOUR1:DATA?'); %read data from buffer
% Convert values to numbers.% First character in string is “{“
% and 2 latest are empty spaces and last is “}”.

signal_num=str2num(signal_str(1,2:length(signal_str)-3)); %str2num par default pas str2double
sample_number = length(signal_num)-1;
x_axis_raw= (0:1:sample_number);
fs=15.258e3; %sampling rate/frequency See decimation factor !!!
x_axis= x_axis_raw./fs; %sample number / Fs, for Fs see decimation factor
%FFT
FFT_signal=fft(signal_num);
FFT_abs=abs(FFT_signal);
N = length(FFT_abs);
fgrid = fs*(0:(N-1))/(N);
FFT_abs = FFT_abs(100:floor(N/4));
fgrid = fgrid(100:floor(N/4));

% % % % % %
% >> Fs=100;  % sample @ 100 Hz
% >> T=10;    % collect data long enough for at least a couple cycles
% >> N=T*Fs;  % compute number samples needed for above
% >> t=linspace(0,10,N);  % and generate time vector
% >> y = 0.7*sin(2*pi*4*t)+randn(size(t));  % and a sample signal around 4Hz
% >> Y = fft(y)/N;              % FFT
% >> PSD=2*abs(Y(1:L/2+1));     % and the PSD one-sided
% >> f=linspace(0,Fs/2,N/2+1);  % compute freq vector for Fs
% >> plot(f,PSD)                % plot the result
% >> hold all                   % we're going to put another on top...
% >> [~,ix]=max(PSD)  % and where's the peak located???
% ix =
%   41
% >> f(ix)
% ans =
%    4
% >> Fs=250;  % ok, let's change sample rate
% >> N=T*Fs;  % next steps just repeat for new sampling vector...
% >> t=linspace(0,10,N);
% >> y = 0.7*sin(2*pi*4*t)+randn(size(t));
% >> Y = fft(y)/N;
% >> PSD=2*abs(Y(1:N/2+1));
% >> f=linspace(0,Fs/2,N/2+1);
% >> plot(f,PSD)
% >> xlim([0 20])  % move scale down to show peak area more clearly
% % % % % %

%plot
figure
plot(x_axis,signal_num)
grid on
ylabel('Voltage(V)')
xlabel('Time(s)')
%hold off
figure
plot(fgrid,FFT_abs);
ylabel('power ($V/\sqrt(Hz)$)','Interpreter','latex')
xlabel('frequencies')
grid on

% }
%sampling rate/data_points
% 
% % Plot results
%     plot_data_1D(HW, data_1D);
% 
%     amplitude2rmsUcoil = data(1).Amplitude2Uin(1) / HW.RX(SeqOut.AQ(1).Device).LNAGain / sqrt(2);
% 
%     hf = figure(9); clf(hf);
%     hax = axes(hf);
%     plot(hax, data(1).time_of_tRep, abs(data(1).data)*amplitude2rmsUcoil*1e6);
%     title(hax, 'Acquired signal');
%     xlabel(hax, 'time in s');
%     ylabel(hax, sprintf('RMS voltage at coil in %cV', char(181)));
%     grid(hax, 'on');
% 
%     % hf = figure(10); clf(hf);
%     % hax = axes(hf);
%     % plot(hax, data.time_of_tRep, angle(data.data));
%     % title(hax, 'Phase of acquired signal');
%     % ylabel(hax, 'phase in rad');
%     % xlabel(hax, 'Time in s');
%     % grid(hax, 'on');
% 
%     hf = figure(11); clf(hf);
%     hax = axes(hf);
%     % without CIC filter correction
%     plot(hax, ...
%       data(1).f_fft1_data, ...
%       abs(data(1).fft1_data)./data(1).cic_corr*amplitude2rmsUcoil*1e6);
%     % % with CIC filter correction
%     % plot(hax, data(1).f_fft1_data, abs(data(1).fft1_data).*amplitude2rmsUcoil*1e6);
%     xlim(hax, [data(1).f_fft1_data(1), data(1).f_fft1_data(end)])
%     title(hax, 'FFT of acquired signal without CIC correction');
%     xlabel(hax, 'frequency in Hz');
%     ylabel(hax, sprintf('RMS voltage at coil in %cV', char(181)));
%     grid(hax, 'on');
% 
%     hf = figure(12); clf(hf);
%     hax = axes(hf);
%     % plot(hax, data(1).f_fft1_data-SeqOut.AQ(1).Frequency(1),abs(data(1).fft1_data)*amplitude2rmsUcoil*1e6); xlabel(hax, 'Frequency in Hz');  % offset frequency
%     plot(hax, (data(1).f_fft1_data/SeqOut.AQ(1).Frequency(1)-1)*1e6, ...
%       abs(data(1).fft1_data)*amplitude2rmsUcoil*1e6);
%     title(hax, 'FFT of acquired signal with CIC correction');
%     xlabel(hax, 'Frequency in ppm'); % offset ppm
%     ylabel(hax, sprintf('RMS voltage at coil in %cV', char(181)));
%     grid(hax, 'on');
    
%% Close connection with Red Pitaya
fclose(tcpipObj);
