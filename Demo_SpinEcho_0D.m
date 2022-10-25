%% Demo Sequence "Simple Spin Echo 0D"
% This demo sequence demonstrates how to create a basic Spin Echo.

%% Simple Spin Echo 0D
% Preparations
    LoadSystem                                      % load system parameters
    [HW,mySave]  = Find_Frequency_Sweep(HW, mySave, 10, [], 1); % find magnet frequency

% Define parameters required for the measurement

    % Parameters used for timing calculations
    Seq.tEcho       =   50e-3;                      % echo time 50e-3
    Seq.p90         =   45e-6; % or HW.tFlip90Def;  % duration of 1st TX pulse 45e-6
    Seq.p180        =   90e-6; % or HW.tFlip180Def; % duration of 2nd TX pulse 90e-6
    Seq.plotSeq     =   [];                         % plot sequence off

    % Sequence parameters
    Seq.tRep        =   100e-3;                     % repetition time


    % RF transmission parameters
    TX.Start        = [ -Seq.p90/2; ...             % start time of 1st TX pulse
                        Seq.tEcho/2-Seq.p180/2];    % Start Time of 2nd TX pulse
    TX.Duration     = [ Seq.p90; ...                % duration of 1st TX pulse
                        Seq.p180];                  % duration of 2nd TX pulse
    TX.Frequency    = [ HW.fLarmor; ...             % frequency of 1st TX pulse
                        HW.fLarmor];                % frequency of 2nd TX pulse
    TX.Phase        = [ 0; ...                      % phase of 1st TX pulse
                        -90];                       % phase of 2nd TX pulse

    % Acquisition parameters
    AQ.Start        = [ 100e-6; ...                 % acquisition start 1st AQ window
                        Seq.tEcho/2+500e-6];        % acquisition start 2nd AQ window
    AQ.fSample      = [ 30e3; ...                   % sampling rate of 1st AQ window
                        50e3];                      % sampling rate of 2nd AQ window
    AQ.nSamples     = [ 128; ...                    % number of samples in 1st AQ window
                        2048];                      % number of samples in 2nd AQ window
    AQ.Frequency    = [ HW.fLarmor; ...             % frequency of 1st AQ window
                        HW.fLarmor];                % frequency of 2nd AQ window
    AQ.Phase        = [ 0; ...                      % phase of 1st AQ window
                        0];                         % phase of 2nd AQ window

% Start measurement
    [ Raw, SeqOut, data, data_1D ] = set_sequence(HW, Seq, AQ, TX, Grad);

% Plot results
    plot_data_1D(HW, data_1D);

%% ------------------------------------------------------------------------
% (C) Copyright 2012-2017 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%------------------------------------------------------------------------
