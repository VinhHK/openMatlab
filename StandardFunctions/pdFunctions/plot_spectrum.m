function spectrumData = plot_spectrum(HW, hf, ppm, spectrumData, waitForClose)
%% Plot spektrum and drag to shift phase
%
%   spectrumData = plot_spectrum(HW, hf, ppm, spectrumData, waitForClose)
%
% INPUT:
%   HW              HW structure
%   hf              Handle to figure with spectrum
%   ppm             x-axis
%   spectrumData    Data with complex spectrum in T
%   waitForClose    If true, blocks execution until figure is closed.
%
% OUTPUT:
%   spectrumData    If waitForClose is true, the modified complex spectrum
%                   in T is returned when the figure is closed. You can
%                   get the current complex spectrum from the open figure
%                   with handle "hf" by the following:
%                       spectrumDataCorr = getappdata(hf, 'spectrumData');
%
% Hold mouse button and move left-right for linear phase, or up-down for
% phase offset.
% Hold "Shift" for faster changes. Hold "Ctrl" for more sensitive changes.
% Double-click to reset phase.
%

%% default input
if nargin < 5
    waitForClose = false;
end

%%
if waitForClose
    waitfor(doPlot(HW, hf, ppm, spectrumData));
else
    doPlot(HW, hf, ppm, spectrumData);
end

%% nested functions to keep spectrumData up-to-date in the main function
    function hf = doPlot(HW, hf, ppm, spectrumData)
        %% Initialize figure with callbacks and initial plot of data
        %
        %   hf = doPlot(HW, hf, ppm, spectrumData)
        %
        % INPUT:
        %   HW              HW structure
        %   hf              Handle to figure with spectrum
        %   ppm             x-axis
        %   spectrumData    Data with complex spectrum
        %
        % OUTPUT:
        %   hf              Handle to figure with spectrum
        %

        %% initialize figure
        hf = figure(hf);
        clf(hf);
        set(hf, 'Name', 'Spectrum - Phase Shift');
        % jhf = get(hf, 'JavaFrame');
        % jj=jhf.getFigurePanelContainer.getComponent(0);
        % jj.setToolTipText('Hold mouse button and move left-right for linear phase, or up-down for phase offset')

        %% save data in figure
        setappdata(hf, 'ppm', ppm);
        setappdata(hf, 'spectrumData', spectrumData);
        setappdata(hf, 'spectrumDataOrig', spectrumData);

        %% plot to axes
        % axes for ticks in Hz (in background)
        haxes(4) = axes('Parent', hf, 'Position', [0.1, 0.1, 0.83, 0.85]);

        % main axes
        % real part of spectrum
        haxes(1) = axes('Parent', hf, 'Position', [0.1, 0.4, 0.83, 0.55]);
%         hReal = plot(haxes(1), ppm, real(spectrumData)/HW.RX.AmplitudeUnitScale);
        hReal = plot(haxes(1), ppm, real(spectrumData));
        hold(haxes(1), 'all');
        setappdata(hf, 'hReal', hReal);
        % absolute
%         hAbs = plot(haxes(1), ppm, abs(spectrumData)/HW.RX.AmplitudeUnitScale, ':', 'HitTest', 'off');
        hAbs = plot(haxes(1), ppm, abs(spectrumData), ':', 'HitTest', 'off');
        lineColor = get(hAbs, 'Color');
        has_hg2 = [100, 1] * sscanf(version, '%d.', 2) >= 804; % Before Matlab 8.4, transparency does not work
        transp = .8;
        if has_hg2
            lineColor(4) = transp;
        else
            lineColor = (transp*lineColor + (1-transp)*[1 1 1]); % Background is white
        end
        set(hAbs, 'Color', lineColor);
        % imaginary part
%         hImag = plot(haxes(1), ppm, imag(spectrumData)/HW.RX.AmplitudeUnitScale, ':', 'HitTest', 'off');
        hImag = plot(haxes(1), ppm, imag(spectrumData), ':', 'HitTest', 'off');
        lineColor = get(hImag, 'Color');
        transp = .8;
        if has_hg2
            lineColor(4) = transp;
        else
            lineColor = (transp*lineColor + (1-transp)*[1 1 1]); % Background is white
        end
        set(hImag, 'Color', lineColor);
        setappdata(hf, 'hImag', hImag);
        % axes
        ylims = get(haxes(1), 'YLim');
%         ylabelStr = [HW.RX.AmplitudeName, ' in ' HW.RX.AmplitudeUnit];
        ylabelStr = ['amplitude', ' relative'];
        ylabel(haxes(1), ylabelStr);
        set(haxes(1), 'XTickLabel', [], 'YLim', [-ylims(2)*.5, ylims(2)]);
        legend({'real', 'abs', 'imag'})

        % axes for integral
        haxes(2) = axes('Parent', hf, 'Position', [0.1, 0.25, 0.83, 0.15]);
        acc_spec = -sign(diff(ppm(1:2)))*(cumsum(real(spectrumData))/sum(real(spectrumData))-.5)+.5;
        hInt = plot(haxes(2), ppm, acc_spec);
        hold(haxes(2), 'all');
        acc_spec((acc_spec>=0)&(acc_spec<=1)) = NaN;
        hIntOff = plot(haxes(2), ppm, acc_spec, 'HitTest', 'off');
        set(haxes(2), 'Color', 'none', 'YLim', [-.1 1.1], 'XTickLabel', []);
        set(hInt, 'Color', [0 0 .5], 'LineWidth', 0.2);
        set(hIntOff, 'Color', 'r', 'LineWidth', 0.2);
        setappdata(hf, 'hInt', hInt);
        setappdata(hf, 'hIntOff', hIntOff);
        ylabel(haxes(2), 'spectrum integral');

        % axes for phase
        haxes(3) = axes('Parent', hf, 'Position', [0.1, 0.1, 0.83, 0.15]);
        hPhase = plot(haxes(3), ppm, angle(spectrumData));
        labelPpm = xlabel(haxes(3), 'ppm');
        set(labelPpm, 'Units', 'normalized');
        posLabelPpm = get(labelPpm, 'Position');
        posLabelPpm(1:2) = [1.025,-0.025];
        set(labelPpm, 'Position', posLabelPpm, 'HorizontalAlignment', 'left');
        ylabel(haxes(3), 'phase in rad');
        setappdata(hf, 'hPhase', hPhase);
        set(haxes(3), 'YLim', [-pi pi]);

        linkaxes(haxes(1:3), 'x');
        set(haxes(1:3), 'XGrid', 'on', 'YGrid', 'on', 'Color', 'none');
        set(haxes, 'XDir', 'reverse');
        setappdata(hf, 'haxes', haxes);
        
        % Add button down functions
        set(haxes(1:3), 'ButtonDownFcn', ...
            @(hObject, eventdata) axes1_ButtonDownFcn(hObject, eventdata, hf, HW));
        
        % settings for axes with Hz
        xlims = get(haxes(1), 'XLim');
        set(haxes(4), 'XLim', xlims/1e6*HW.fLarmor, 'YTick', [], 'XGrid', 'on', ...
            'GridLineStyle', ':', 'XAxisLocation', 'top', 'HitTest', 'off');
        transp = .9;
        gridColor = [.2 .8 .3];
        if has_hg2
            set(haxes(4), 'GridColor', gridColor, 'GridAlpha', transp);
        else
            gridColor = (transp*gridColor + (1-transp)*[1 1 1]); % Background is white
            set(haxes(4), 'XColor', gridColor, 'Box', 'on');
        end
        labelHz = xlabel(haxes(4), 'Hz');
        set(labelHz, 'Units', 'normalized');
        posLabelHz = get(labelHz, 'Position');
        posLabelHz(1:2) = [1.025,1];
        set(labelHz, 'Position', posLabelHz, 'HorizontalAlignment', 'left');
       
        hZoom = zoom(hf);
        set(hZoom, 'ActionPostCallback', @(hObject, event) spectrum_ActionPostCallback(hObject, event, haxes, HW));
        hPan = pan(hf);
        set(hPan, 'ActionPostCallback', @(hObject, event) spectrum_ActionPostCallback(hObject, event, haxes, HW));
    end

    function spectrum_ActionPostCallback(hObject, eventdata, haxes, HW)
        xlims = get(haxes(1), 'XLim');
        set(haxes(4), 'XLim', xlims/1e6*HW.fLarmor);
    end

    function axes1_ButtonDownFcn(hObject, eventdata, hf, HW)
        %% Executes on Mouse button click down
        %
        %   axes1_ButtonDownFcn(hObject, eventdata, hf)
        %
        % INPUT:
        %   hObject     handle to axes1
        %   eventdata   contains info to mouse action
        %   hf          handle to figure
        %
        % Sets WindowButtonMotionFcn and WindowButtonUpFcn for the figure.
        % Resets phase on double-click.

        if strcmpi(get(hf, 'SelectionType'), 'open') % double-click
            % reset data
            spectrumData = getappdata(hf, 'spectrumDataOrig');
            setappdata(hf, 'spectrumData', spectrumData);
            updatePlot(hf, ppm, spectrumData, HW)
        else
            if strcmpi(get(hf, 'SelectionType'), 'extend') % shift+click
                sensitivity = .2;
            elseif strcmpi(get(hf, 'SelectionType'), 'alt') % ctrl+click
                sensitivity = 5;
            else % normal click
                sensitivity = 1;
            end
            apos = get(hObject, 'CurrentPoint'); % in axes coordinates
            clickPosition = apos(1,1:2);
            set(hf, 'WindowButtonMotionFcn', @(hObject, eventdata) spectrum_WindowButtonMotionFcn(hObject, eventdata, clickPosition, sensitivity, HW));
            set(hf, 'WindowButtonUpFcn', @(hObject, eventdata) spectrum_WindowButtonUpFcn(hObject, eventdata));
        end
    end

    function spectrum_WindowButtonUpFcn(hObject, eventdata)
        %% Executes on mouse button release
        %
        %   spectrum_WindowButtonUpFcn(hObject, eventdata)
        %
        % INPUT:
        %   hObject     handle to figure
        %   eventdata   contains info to mouse action
        %
        % Clears WindowButtonMotionFcn and WindowButtonUpFcn for the figure.

        set(hObject, 'WindowButtonMotionFcn', []);
        set(hObject, 'WindowButtonUpFcn', []);
        spectrum_WindowButtonMotionFcn([]);
    end

    function spectrum_WindowButtonMotionFcn(hObject, eventdata, clickPosition, sensitivity, HW)
        %% Executes on mouse pointer movement
        %
        %   spectrum_WindowButtonMotionFcn(hObject, eventdata, clickPosition)
        %
        % INPUT:
        %   hObject         handle to figure
        %   eventdata       contains info to mouse action
        %   clickPosition   coordinates of mouse button down in axes coordinates
        %
        % Actual correction of phase on mouse movement.
        % Initialize by calling function with one empty argument.

        %% initialize
        persistent oldpos

        %% early return to initialize oldpos
        if isempty(hObject)
            oldpos = [];
            return
        end

        %% mouse action
        apos = get(hObject, 'CurrentPoint'); % in pixels
        if isempty(oldpos)
            oldpos = apos;
            return
        end
        pixelMoved = apos-oldpos;
        oldpos = apos;

        spectrumData = getappdata(hObject, 'spectrumData');
        % constant phase
        spectrumData = spectrumData*exp(-1i*pixelMoved(2)/4000/sensitivity*2*pi);
        % linear phase
        ppm = getappdata(hObject, 'ppm');
        linearPhase = pixelMoved(1)/8000/sensitivity*2*pi * (ppm - clickPosition(1));
        spectrumData = spectrumData.*exp(-1i*linearPhase);
        % update data
        setappdata(hObject, 'spectrumData', spectrumData);
        updatePlot(hObject, ppm, spectrumData, HW)
    end
end

function updatePlot(hObject, ppm, spectrumData, HW)
%% update lines in plots with spectrumData

hReal = getappdata(hObject, 'hReal');
% set(hReal, 'YData', real(spectrumData)/HW.RX.AmplitudeUnitScale);
set(hReal, 'YData', real(spectrumData));
hImag = getappdata(hObject, 'hImag');
% set(hImag, 'YData', imag(spectrumData)/HW.RX.AmplitudeUnitScale);
set(hImag, 'YData', imag(spectrumData));
acc_spec = -sign(diff(ppm(1:2)))*(cumsum(real(spectrumData))/sum(real(spectrumData))-.5)+.5;
hInt = getappdata(hObject, 'hInt');
set(hInt, 'YData', acc_spec);
acc_spec((acc_spec>=0)&(acc_spec<=1)) = NaN;
hIntOff = getappdata(hObject, 'hIntOff');
set(hIntOff, 'YData', acc_spec);
hPhase = getappdata(hObject, 'hPhase');
set(hPhase, 'YData', angle(spectrumData));

end