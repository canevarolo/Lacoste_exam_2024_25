% A photovoltaic panel receives an average daily irradiance power of 500 W/m² in the spring season.
% Considering the air, with which the panel exchanges heat on its rear side with a heat transfer
% coefficient of 10 W/m²K, at a temperature of 20 °C, calculate the temperature distribution along
% the thickness of the panel, neglecting other components. In the calculation, consider the thickness
% of the protective glass, with a thermal conductivity of 0.76 W/mK and a thickness of 3.5 mm,
% while the silicon has a thickness of 0.5 mm and a thermal conductivity of 148 W/mK.

% Laib 4, exercise 4
% Simone Canevarolo
% S269893
% 6/11/2024

clear all
close all
clc

qirr = 500; % irradiance flow, W/m^2
henv = 10; % W/m/K, environment
Tenv = 293; % K, environment

kkglass = 0.76; % W/m/K
ssglass = 3.5e-3; % m

kksil = 148; % W/m/K
sssil = 0.5e-3; % m

dx = 1e-5;
xx = (0:dx:(ssglass+sssil))';
Nx = length(xx);

xxglass = (0:dx:ssglass)';
Nxglass = length(xxglass);

sub_diag = ones(Nx,1);
main_diag = -2*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = zeros(Nx,1);

AA(1,1) = 1;
AA(1,2) = -1;
bb(1) = qirr*dx/kkglass;

AA(Nxglass,Nxglass-1) = kkglass;
AA(Nxglass,Nxglass) = -kksil-kkglass;
AA(Nxglass,Nxglass+1) = kksil;

AA(end,end-1) = -1;
AA(end,end) = 1+henv*dx/kksil;
bb(end) = Tenv*henv*dx/kksil;

TT = AA\bb;

figure(1)
plot(xx,TT-273,'r','LineWidth',2)
title('Temperature along the panel')
xlabel('Thickness [m]')
ylabel('Temperature [°C]')
legend('Temperature','Location','best')