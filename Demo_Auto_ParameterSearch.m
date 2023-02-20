%% Demo Sequence "Auto Parameter Search"
% automatically find the correct frequency, RF 90 degree pulse duration and
% shim values; values will be appended to "User\MagnetShimCal.m" and
% "PaUout2AmplitudeCal.m"
% Use 10 mm Oil sample for best results

%% Auto Parameter Search
% Preparations
LoadSystem;                           % load system parameters

%% Find_Frequency_Sweep( HW, mySave, maxtime, span, doplot, tPulse90)
% parameters:
%  HW, mySave; these parameters are created automatically by LoadSystem
LoadSystem;                           % load system parameters

minTime               = 1;            % minimum time in seconds since the last Find_Frequency; i.e. if set to 10, find freq will only be run, if 10 seconds have passed since the last search
span                  = 400e3;        % searching span of frequency around fLarmor (HW.B0*HW.FindFrequencyGamma/2/pi) 400e3
doplot                = 1;            % plot sequence and data
tPulse90              = [];           % duration of the 90 degrees pulse ([] uses HW.FindFrequencySweep.tPulse90 or HW.tFlip90)
nMeasurements         = 21;           % number of measurements (TX and AQ trains)
nSamples              = 100;          % number of samples in each measurement

[HW, mySave] = Find_Frequency_Sweep(HW, mySave, minTime, span, doplot, tPulse90, nMeasurements, nSamples);

%% Find_Shim(HW, mySave, maxtime, doplot, iterations, tEcho, T1, ShimStart, ShimStep)
% parameters:
%  HW, mySave; these parameters are created automatically by LoadSystem
LoadSystem;                           % load system parameters
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);

minTime               = 10;           % minimum time in seconds since the last Find_Frequency_Sweep.10
doplot                = 1;            % plot sequence and data 1
Seq.iterations        = 100;          % number of iterations used for shim search
Seq.tEcho             = 100e-3;        % echo time 50e-3
Seq.RepetitionTime    = 1;          % repetition time: ~estimated T1 value of the used sample * 3 e.g. 0.5 s (Seq.RepetitionTime=3*T1)
Seq.ShimStart         = [];           % start with these 1x4 shim values ([] for defaults)
Seq.ShimStep          = [];           % 1x4 vector with step widths ([] for defaults)
Seq.nEchos            = 1;            % use FID or Echo for shimming (0 for Fid 1,2,3... for nEcho)
Seq.use_nEchos_Frequency = 1;         % find frequency of nEchos

Seq.useSliceSelect    = 0;            % use slice gradient
clear SliceSelect
SliceSelect.MaxGradAmpSlice = 0.05;   % maximum gradient amplitude during slice excitation
SliceSelect.alfa      = 0.0*pi;       % rotation around x axis
SliceSelect.phi       = 0.0*pi;       % rotation around y axis
SliceSelect.theta     = 0.5*pi;       % rotation around z axis
SliceSelect.thickness = 0.008;        % slice thickness

[HW, mySave] = Find_Shim(HW, mySave, minTime, doplot, Seq, SliceSelect);

LoadSystem;                                 % load new parameters from file.
disp(['New shim values: ' num2str(HW.MagnetShim) ' T/m.']);

% search frequency with new pulse value
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);

%% Find_PulseDuration(HW, mySave, maxtime, doplot, iterations, tPuls90, T1)
% parameters:
%  HW, mySave; these parameters are created automatically by LoadSystem
LoadSystem;                           % load system parameters

minTime               = 0;            % minimum time in seconds since the last Find_Pulseduration.
doplot                = 1;            % plot sequence and data
iterations            = 1;            % number of iterations used for pulse length determination / 0: Set the tPulse90 without searching.
tPulse90              = [];           % duration of the 90 degrees pulse
T1                    = 0.15;         % estimated T1 value of the used sample 0.15

[HW, mySave] = Find_PulseDuration(HW, mySave, minTime, doplot, iterations, tPulse90, T1);

LoadSystem;                           % load new parameters from file.
disp(['New pulse duration: ' num2str(HW.tFlip90Def) ' s.']);


%% Find_Frequency_FID(HW, mySave, minTime, iterations, tAQ)
% parameters:
%  HW, mySave; these parameters are created automatically by LoadSystem
LoadSystem;                           % load system parameters
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 100);

minTime               = 0;            % minimum time in seconds since the last Find_Frequency.
iterations            = 2;            % number of iterations used for freqeuency determination
tAQ                   = 2e-3;         % acquisition time in s
tAQStart              = [];           % start time of acquisition window in s
oldDoPlot = HW.FindFrequencyFID.doPlot;
HW.FindFrequencyFID.doPlot = 1;       % plot FID

[HW, mySave] = Find_Frequency_FID(HW, mySave, minTime, iterations, tAQ, tAQStart);

HW.FindFrequencyFID.doPlot = oldDoPlot;


%% ----------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%------------------------------------------------------------------------
