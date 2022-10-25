f=25e6;
w0=2*pi*f;
dcoil=1e-3;

% Wire_dia=50e-6:10e-6:500e-6;
% nturn=1:1:50;

% lcoil=nturn.*Wire_dia;

% Bxy=(nturn.*4.*pi.*1e-7)./(dcoil.*(sqrt(1+(lcoil/dcoil).^2)));
% Signal=(w0^2.*Bxy)./sqrt(2);


for i = 1:length(Wire_dia)
    for j=1:length(nturn)
        nturn=nturn(j);
        Bxy=(nturn(j).*4.*pi.*1e-7)./(dcoil.*(sqrt(1+(lcoil/dcoil).^2)));
        [Data(i,j)]=(w0^2.*Bxy)./sqrt(2);
    end
end
% Pour récupérer les valeurs à plotter
for i=1:length(Data)
    for j=1:length(Data)
        X(i,j) = Data(i,j).X;
        Y(i,j) = Data(i,j).Y;
        Z(i,j) = Data(i,j).Z;
    end
end

figure;

surf(X,Y,Z);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Qz=f(\alpha,n)');
grid on;
% 
% figure
% plot(lcoil*1e3,Signal);
% ylabel('Sensitivity Bxy ($T$)','Interpreter','latex')
% xlabel('Coil length (mm)')
% zlabel('Number of turns')
% grid on
