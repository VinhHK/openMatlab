function hf = plot_Reflection(HW, Network, hax)
%% Plot reflection S11 (amplitude and phase) in dB of network
%
%     hf = plot_Reflection(HW, Network, hax)
%
% Annotations are added to the graphs marking the Larmor frequency and the
% resonance frequency of the LC circuit. If possible, the Q value is calculated
% using the -3 dB bandwidth and shown in the annotation.
%
% INPUT:
%   HW
%           HW object or structure.
%   Network
%           Network structure as returned by sequence_Network.
%   hax
%           Optional handle to axes that is used for plotting the reflection
%           S11 in dB. If omitted, a figure is created that displays two axes
%           with the amplitude in dB and the phase of the reflection S11.
%
% OUTPUT:
%   hf
%           Only if called with two input argument or if hax is empty, the
%           handle to the figure that contains the reflection graphs.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% default input
if nargin == 2, hax = []; end

%% create figure
if isempty(hax)
  hf = figure(200); clf(hf);
  set(hf, 'Name', 'Reflection and Phase', 'NumberTitle', 'off');
  hax = subplot(2,1,1, 'Parent', hf);
end

%% plot reflection S11
% plot(hax, Network.Frequency,20*log10(abs(Network.Reflection)), 'LineWidth' ,2)
% ylim(hax, [floor((min(20*log10(abs(Network.Reflection)))-2)/5)*5, ceil(max(20*log10(abs(Network.Reflection)))/5)*5])
Frequency = interpn(Network.Frequency, 4);
Reflection = interpn(Network.Reflection, 4, 'spline');
[MinReflection, iMinReflection] = min(double(Reflection));
fMinReflection = Frequency(iMinReflection);
plot(hax, Frequency, 20*log10(abs(Reflection)), 'LineWidth', 2);
ydiff = diff([floor((min(20*log10(abs(Reflection))))/5)*5, ceil(max(20*log10(abs(Reflection)))/5)*5]);
ylim(hax, [floor((min(20*log10(abs(Reflection)))-ydiff/10)/5)*5, ceil(max(20*log10(abs(Reflection)))/5)*5]);
xlim(hax, [min(Frequency), max(Frequency)]);
ylabel(hax, 'Reflection S11 in dB');
xlabel(hax, 'Frequency in Hz');
grid(hax, 'on');

% add labels to the graph
% Larmor frequency
reflectionAtFLarmor = double(20*log10(abs(interp1(Frequency, Reflection, HW.fLarmor, 'nearest'))));
if fMinReflection < HW.fLarmor
  text(HW.fLarmor, ...
    reflectionAtFLarmor, ...
    '\uparrow', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'Rotation', 0, 'Parent', hax);
    % '\rightarrow', ...
    % 'HorizontalAlignment', 'right', ...
    % 'Rotation', 90, 'Parent', hax, ...
    % 'Color', [0.106 0.31 0.208]);

  text(HW.fLarmor, ...
    reflectionAtFLarmor, ...
    sprintf('   %.1f dB @ f_L', reflectionAtFLarmor), ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', ...
    'Rotation', 0, 'Parent', hax, ...
    'Color', [0.106, 0.31, 0.208])
else
  text(HW.fLarmor, ...
    reflectionAtFLarmor, ...
    '\uparrow', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'Rotation', 0, 'Parent', hax);
    % '\leftarrow', ...
    % 'HorizontalAlignment', 'left', ...
    % 'Rotation', 90, 'Parent', hax, ...
    % 'Color', [0.106 0.31 0.208]);

  text(HW.fLarmor, ...
    reflectionAtFLarmor, ...
    sprintf('%.1f dB @ f_L ', reflectionAtFLarmor), ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', ...
    'Rotation', 0, 'Parent', hax, ...
    'Color', [0.106, 0.31, 0.208]);
end

text(fMinReflection, ...
  20*log10(abs(MinReflection)), ...
  '\uparrow', ...
  'FontWeight', 'bold', ...
  'HorizontalAlignment', 'center', ...
  'VerticalAlignment', 'top', ...
  'Rotation', 0, 'Parent', hax);

% Q value
Qstr = '';
if min(20*log10(abs(Reflection))) < -4
  if max(20*log10(abs(Reflection(1)))) > -2.9
    if max(20*log10(abs(Reflection(end)))) > -2.9
      StartIndex = find(20*log10(abs(Reflection))<-3, 1, 'first');
      StopIndex = find(20*log10(abs(Reflection))<-3, 1, 'last');
      [~, MinIndex] = min(20*log10(abs(Reflection)));
      Qstr = sprintf(' Q = %.1f', Frequency(MinIndex)./(Frequency(StopIndex)-Frequency(StartIndex)));
      % text(fMinReflection, ...
      %   20*log10(abs(MinReflection)), ...
      %   sprintf('Q = %.1f \\rightarrow', Frequency(MinIndex)./(Frequency(StopIndex)-Frequency(StartIndex))), ...
      %   'HorizontalAlignment', 'right', ...
      %   'VerticalAlignment', 'middle', ...
      %   'Rotation', 0, 'Parent', ah)
    end
  end
end

text(fMinReflection, ...
  20*log10(abs(MinReflection)), ...
  sprintf(['\n%.1f dB @ %.3f MHz' Qstr], 20*log10(abs(MinReflection)), fMinReflection/1e6), ...
  'HorizontalAlignment', 'center', ...
  'VerticalAlignment', 'top', ...
  'FontWeight', 'bold', ...
  'Rotation', 0, 'Parent', hax);

if exist('hf', 'var')
  % phase of reflection
  hax2 = subplot(2,1,2, 'Parent', hf);
  plot(hax2, Frequency, angle(Reflection), 'LineWidth', 2);
  xlim(hax2, [min(Frequency), max(Frequency)]);
  ylabel(hax2, {'Angle', 'rad'});
  xlabel(hax2, 'Frequency Hz');
  text(HW.fLarmor, ...
    angle(interp1(Frequency,Reflection, HW.fLarmor, 'nearest')), ...
    '\leftarrow fLarmor', ...
    'HorizontalAlignment', 'left', ...
    'Rotation', 90, 'Parent', hax2);
end

% drawnow
drawnow('expose');

end
