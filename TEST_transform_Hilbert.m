clc
close all
clear all

open(['Test_Releve_NMR/FID_TEST.fig']);
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
zdata = get(a, 'ZData');

%% values to be use without Hilbert IQ Transform
f_LO=24374235+2000;
Fs=15.258e3;    % Sampling frequency
T = 1/Fs;   % Sampling period
L = length(ydata);             % Length of signal
t = xdata;        % Time vector
fgrid=Fs/L*(-L/2:L/2-1);
fid = ydata; %;sin(2*pi*30*t);

%% Hilbert transform

fid_complex = hilbert(fid);
% fid_complex = smoothdata(fid_complex); %auto smoothing
%fgrid = ((fgrid*1e6)/(f_LO))-60; %Hz to PPM with recalibration of the 0ppm

%% Plotting time
plot(t,real(fid_complex));
hold on
plot(t,imag(fid_complex));
plot(t,abs(fid_complex));
hold off
legend('Real Part','Imaginary Part','Module Part');

%% FFT
Y_real = real(fftshift(fft((fid_complex))));
Y_imag = imag(fftshift(fft((fid_complex))));
Y_abs = abs(fftshift(fft((fid_complex))));

% P2 = (Y);
%P3 = (Y.*exp(-j.*2.*pi.*fgrid.*(4./fgrid)));
% P3 = (Y.*exp(-j.*2.*pi.*fgrid.*(L./Fs)));
% Pz = complex(Y_real+Y_imag);

%% Plot the FFT
% figure
% 
%  subplot(2,1,1);
% plot(fgrid,(Y_real),fgrid,imag(Y_real),fgrid,((Y_abs)))
% grid on
% legend('Real Hilbert','imag Hilbert','abs Hilbert')
% ylabel('Voltage(V)')
% xlabel('Freq(Hz)')
% 
%  subplot(2,1,2);
% plot(fgrid,imag(Y_imag))
% grid on
% title('Imag Hilbert')
% ylabel('Voltage(V)')
% xlabel('Freq(Hz)')

f_Tx = 24372832;
f_LO = f_Tx+2000; 
fgrid = ((fgrid*1e6)/(f_LO))-60;
t=0:0.0001:1;
figure

plot(fgrid,(Y_abs))
grid on
legend('FFT')
ylabel('Voltage(V)')
xlabel('Chemical Shift (PPM)')

%% LOD selon Trejo Rosillo
%Sensibilite massique (Kg/mol) = RSB/Concentration massique (mol.kg-1)
%Sensibilite molaire (mol-1) = RSB/Concentration molaire (nombre de mol)
% nombre de mole = masse(gramme)/masse molaire (gramme par mole)
% LoD = 3/Sensibilité  (dépend de la sensibilité utilisé)

Somme_bruit=0; %integrale sur 14000 a 16000
for i=14000:16000
    Somme_bruit=Somme_bruit+Y_abs(i);
end

Somme_bruit=Somme_bruit/(16000-14000); %moyenne bruit

Somme_bruit= Somme_bruit*(max(fgrid)-min(fgrid)); %integrale du bruit total

%[pks,locs,widths,proms] = findpeaks(Y_abs,fgrid);

[pks,locs,widths,proms] = findpeaks(Y_abs,fgrid,'MinPeakProminence',4,'Annotate','extents');
% Modelisation en lorentzienne du signal 
% findpeaks(Y_abs,fgrid,'MinPeakProminence',4,'Annotate','extents');

%% Lorentzienne
T2=sqrt(3)./(pi*widths);

Lz_eq= ((pks)/(2.*pi.*T2))./sqrt(((1./(2*pi*T2)).^2) + (fgrid-locs).^2 );

hold on
plot(fgrid,Lz_eq)

Somme_Lz=0;
for i=1:length(fgrid)
    Somme_Lz=Somme_Lz+Lz_eq(i);
end

RSB = Somme_Lz/Somme_bruit;

%Concentration 
%Mineral oil 0.86g/mL inside a 10mm cylinder i.e. 0.79mL with  452g/mol 
Concentration = 0.68/452; %nombre de mole = gramme / gramme par mole
Sensibilite_molaire = RSB/Concentration;

%concentration volumique en g/mol ??

LoD = 3/Sensibilite_molaire; %LOD en mol
fprintf('Limit of Detection = %d \n',LoD)

%% Figures 3D
% figure
% plot3(fgrid,Y_real,Y_imag)
% 
% figure
% plot3(fgrid,real(fid_complex),imag(fid_complex))
