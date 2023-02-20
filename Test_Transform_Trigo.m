clc
clear all
close all


%x = -10*pi:0.005:10*pi; 

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;     % Length of signal
t = (0:L-1)*T;        % Time vector

y=cos(2*pi*40*t).*exp(-t/0.1);

y2=hilbert(y);

plot(t,imag(y2),t,real(y2),t,abs(y2))
legend('imag Hilbert y','real Hilbert y','abs Hilbert y')
grid on

%% FFT
Y2=fftshift(fft(y2));

Y_real = fftshift(fft(real(Y2)));
Y_imag = fftshift(fft(imag(Y2)));
Y_abs = fftshift(fft(abs(Y2)));

figure
plot(t,Y2)
legend('FFT de Hilbert y')
