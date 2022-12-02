%% Demo Sequence "FID - Free Induction Decay"
% This demo sequence demonstrates how to create an FID.
% It applies an excitation pulse and acquires the FID.

%% FID
% Preparations
    LoadSystem;                                     % load system parameters

    %[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10, [], 1);  % find magnet frequency
    %HW.fLarmor = 24371171.6266625;                              % override Larmor frequency

    % Parameters used for timing calculations
    Seq.p90         = 4e-6;%HW.tFlip90Def                % duration of TX pulse
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

% Start measurement
    [Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);

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
    plot(hax, (data(1).f_fft1_data/SeqOut.AQ(1).Frequency(1)-1)*1e6-233.4, ...
      abs((data(1).fft1_data))*amplitude2rmsUcoil*1e6);
    title(hax, 'FFT of acquired signal with CIC correction');
    xlabel(hax, 'Frequency in ppm'); % offset ppm
    ylabel(hax, sprintf('RMS voltage at coil in %cV', char(181)));
    grid(hax, 'on');

    % select figure with FID
    figure(9);

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
