% open('FID_TEST.fig');
% a = get(gca,'Children');
% xdata = get(a, 'XData');
% ydata = get(a, 'YData');
% zdata = get(a, 'ZData');
clc
close all


[figures, path] = uigetfile('Multiple_acquisition_save/', '*.fig', 'Multiselect', 'on');   
for n = 1:length(figures)
    Multi_Figs = [path, filesep,figures{n}];
    Op = openfig(Multi_Figs,'invisible');
    h = findobj(gca, 'Type', 'line');
    S.data_x(:,n) = get(h, 'Xdata');
    S.data_y(:,n) = get(h,'Ydata');
end

% S = dir('Multiple_acquisition_save/*.txt');
% for k = 1:numel(S)
%     S(k).data = readmatrix(S(k).name);
% end


%%
%f_LO=24372832+2000;
Fs=15.258e3;    % Sampling frequency
T = 1/Fs;   % Sampling period
L = length(S.data_y(n,:)); %L = length(S(1).data(:,2));
fgrid=Fs/L*(-L/2:L/2-1);
%3 Time / 2 Data FID
% plot(PT_FID_0min_fig(:,3),PT_FID_0min_fig(:,2))

%%

% t21 = finddelay(PT_FID_0min_fig(:,3),PT_FID_100min_fig(:,3))
%
% e1 = S(1).data(:,2);
% e2 = S(2).data(:,2);

i=0;
for i=1:n
    if i==1
           y1 = S.data_y(:,1); %S(1).data(:,2);
    end

    [y1,y2] = alignsignals(y1,S.data_y(:,i),Method="xcorr");
    %         remettre dans une variable
    %         resynchro avec une nouvelle
%     y_base = y1;
%     y1=y1(1:15000);
%     y2=y2(1:15000);
%     y1=(y1+y2)/2;
end

%%
figure

for i=1:n
[y1,y2] = alignsignals(y1,S.data_y(:,i),Method="xcorr");
% [y3,y4] = xcorr(e1,e2);

% subplot(2,2,1)
%y2=y2(1900:14000);
plot(y2)
hold on
lgd{i} = strcat('File =',num2str(i)) ;
end
 plot(y1)
 hold off
legend(lgd(1:n));
xlabel("Samples")
ylabel("Amplitude")
hold off
%
%
% subplot(2,2,2)
% plot(y3)
% xlabel("Samples")
% ylabel("Amplitude")
% hold on
% plot(y4)
% hold off
%

%% Hilbert transform

f_LO=24372832+2000;
Fs=15.258e3;    % Sampling frequency
T = 1/Fs;   % Sampling period
L = length(y1); % Length of signal
t = S.data_x(:,1);  % Time vector
fgrid=Fs/L*(-L/2:L/2-1);
fid = y1; %;sin(2*pi*30*t);
fid_complex = hilbert(fid);
%% 


P0 = -1;
P1 = -0.0016;
fgrid=fgrid';

P0_step = 1;
P1_step = 0.001;

%%
figure
% for i=1:1
i=1;
% P1=P1+P1_step;

Y_real = real(fftshift(fft((fid_complex))).*exp(1i*P0).*exp(1i*P1*2*pi*fgrid));

%   Y_imag = imag(fftshift(fft((fid_complex))));
%   .*exp(j+2*pi*fgrid*P1);

plot(fgrid,Y_real);

lgd{i} = strcat('P0=',num2str(P0),'and P1=',num2str(P1)) ;

hold on

% end
f_Tx = 24372832;
f_LO = f_Tx+2000;
fgrid = ((fgrid*1e6)/(f_LO))-60;

grid on
legend(lgd(1:length(lgd)));
ylabel('Voltage(V)')
xlabel('Chemical Shift (PPM)')
hold off