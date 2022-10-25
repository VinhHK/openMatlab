close all
clc
clear all

Id_Ref_Freq=3848;

S_Hi = sparameters('Hi_Z_NDV.s2p')
S_50 = sparameters('50_Ohm_NDV.s2p')
%disp(S)
%rfplot(S)
% HIGH IMPEDANCE %
s11_Hi = rfparam(S_Hi,1,1);
s12_Hi = rfparam(S_Hi,1,2);
s21_Hi = rfparam(S_Hi,2,1);
s22_Hi = rfparam(S_Hi,2,2);
smithplot(S_Hi)
% 50 Ohm %
s11_50 = rfparam(S_50,1,1);
s12_50 = rfparam(S_50,1,2);
s21_50 = rfparam(S_50,2,1);
s22_50 = rfparam(S_50,2,2);
smithplot(S_50)

s12x_Hi = (s12_Hi(Id_Ref_Freq));%@100MHz = line Id_Ref_Freq exact 1.000210880452221e+08 MHz
s22x_Hi = (s22_Hi(Id_Ref_Freq));
s21x_Hi = (s21_Hi(Id_Ref_Freq));
s11x_Hi = (s11_Hi(Id_Ref_Freq));
zo = 50;

s12x_50 = (s12_50(Id_Ref_Freq));%@100MHz = line Id_Ref_Freq exact 1.000210880452221e+08 MHz
s22x_50 = (s22_50(Id_Ref_Freq));
s21x_50 = (s21_50(Id_Ref_Freq));
s11x_50 = (s11_50(Id_Ref_Freq));

Z_Hi = zparameters('Hi_Z_NDV.s2p');

z11_Hi = rfparam(Z_Hi,1,1);
z12_Hi = rfparam(Z_Hi,1,2);
z21_Hi = rfparam(Z_Hi,2,1);
z22_Hi = rfparam(Z_Hi,2,2);

% figure('Name','Z11','NumberTitle','off');
% plot(Z.Frequencies, imag(z11)) ;% this is for imaginary part of Z-, change as per your requrements
% figure('Name','Z12','NumberTitle','off');
% plot(Z.Frequencies, imag(z12)) ;
% figure('Name','Z21','NumberTitle','off');
% plot(Z.Frequencies, imag(z21)) ;
% figure('Name','Z22','NumberTitle','off');
% plot(Z.Frequencies, imag(z22)) ;

 for c= 1:4000
 k(c,:)=sqrt((imag(z12_Hi(c,:)).*imag(z21_Hi(c,:)))/(imag(z11_Hi(c,:)).*imag(z22_Hi(c,:))));
 end
% figure('Name','k factor','NumberTitle','off');
% plot(k)

% Comparaisons des parametres Z obtenus par Matlab et par calcul
Z12_Extracted = (z12_Hi(Id_Ref_Freq)); %@100MHz = line 3750
Z12_Calculated = zo*((2*s12x_Hi)/((1-s11x_Hi)*(1-s22x_Hi)-(s12x_Hi*s21x_Hi)));
Z22_Extracted = (z22_Hi(Id_Ref_Freq)); 
Z22_Calculated = zo*(((1-s11x_Hi)*(1+s22x_Hi)+(s12x_Hi*s21x_Hi))/((1-s11x_Hi)*(1-s22x_Hi)-(s12x_Hi*s21x_Hi)));
Z21_Extracted = (z21_Hi(Id_Ref_Freq)); 
Z21_Calculated = zo*((2*s21x_Hi)/((1-s11x_Hi)*(1-s22x_Hi)-(s12x_Hi*s21x_Hi)));
Z11_Extracted = (z11_Hi(Id_Ref_Freq)); 
Z11_Calculated = zo*(((1+s11x_Hi)*(1-s22x_Hi)+(s12x_Hi*s21x_Hi))/((1-s11x_Hi)*(1-s22x_Hi)-(s12x_Hi*s21x_Hi)));


% Calcul de Gamma %
coefficient_Hi = gammams(S_Hi);
Gamma_Coeff_Hi=real(coefficient_Hi(Id_Ref_Freq));

coefficient_50 = gammams(S_50);
Gamma_Coeff_50 = real(coefficient_50(Id_Ref_Freq));
%Gamma_Hiz= Z1

%Calcul du Z1%
Z1_Hi=zo*((1+Gamma_Coeff_Hi)/(1-Gamma_Coeff_Hi))
Z1_50=zo*((1+Gamma_Coeff_50)/(1-Gamma_Coeff_50))

%Calcul de I%
P=100e-6;
V_RF=sqrt(8*zo*P);
I1_Hi=V_RF/(zo+Z1_Hi);
I1_50=V_RF/(zo+Z1_50);

I_Hi_sur_50=I1_Hi/I1_50

Vfid_Hi_50 = (zo+Z1_50) / (zo+Z1_Hi)



% Calcul du facteur M %

f0 = 100e6; 
wo = 2*pi*f0;
k_100M = k(Id_Ref_Freq); %facteur k à 100MHz
L1 = 78e-9; 
L2 = 78e-9;
R15= 5100; %resistance serie du primaire
R2 = 0.7 ; % Resistance au secondaire

M = imag(z21_Hi(Id_Ref_Freq))/(2*pi*f0)
%M = k_100M*(sqrt(L1*L2))

z2 = R2 + ((M^2*wo^2)/(zo+R15));
z1hi = R15+((M^2*wo^2)/R2)
z150 = R15+((M^2*wo^2)/(2*R2 + ((M^2*wo^2)/(zo+R15))))

Rp_Vfid_Hi_50 = (zo+z150) / (zo+z1hi)


% Calcul des facteurs de qualités %
Q2_P = (L2*wo) / (R2 + ((M^2*wo^2)/(zo+R15)))
Q2 = (L2*wo)/ R2

Vout_hi_50 = Q2 / ((1/2)*sqrt(zo/R2))

% 20mVpp ohm 50
% 2Vpp  Hiz

Vfid_Hi_50_Calculated = (Q2_P*(zo+R15))/((1/2)*sqrt(zo/(R2 + ((M^2*wo^2)/(zo+R15))))*(zo+z150))