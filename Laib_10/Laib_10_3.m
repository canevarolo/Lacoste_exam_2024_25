% A long bar with a thermal conductivity of 1.5 W/mK and a rectangular cross-section of 0.4 m × 0.6 m
% is subject to the following boundary conditions:
% Two faces of the bar are maintained at a uniform temperature of 200°C;
% one face is adiabatic,
% and the remaining face is subject to convective heat transfer with a heat transfer coefficient
% h = 50 W/mK to a fluid at T = 30°C.
% Determine the temperature distribution within the bar and the heat power dissipated to the fluid per unit length.

% Laib 10, exercise 3
% Simone Canevarolo
% S269893
% 08/01/2025

clear all
close all
clc

base = 0.4; % m
alt = 0.6; % m

kk = 1.5; % W/m/K
hh = 50; % W/m^2/K
Tfix = 200+273; % K
Tair = 30+273; % K

dx = 1e-2; % m
dy = dx;

xx = (0:dx:alt)';
yy = (0:dy:base)';

Nx = length(xx);
Ny = length(yy);

Ntot = Nx*Ny;

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

[xmat,ymat] = meshgrid(xx,yy);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

bb = zeros(Ntot,1);


% North

kk = 1:Nx:Ntot-Nx+1;

AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tfix;


% East

jj = Ny;
for ii = 2:Nx-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2-hh/dy;
    AA(kk,kk+1) = 1/2/dx^2;
    bb(kk) = -Tair*hh/dy;

end


% South

kk = Nx:Nx:Ntot;

AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tfix;


% West

jj = 1;
for ii = 2:Nx-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-1) = 1/2/dx^2;
    AA(kk,kk) = -1/dx^2-1/dy^2;
    AA(kk,kk+1) = 1/2/dx^2;
    AA(kk,kk+Nx) = 1/dy^2;

end


TT = AA\bb;

TTT = reshape(TT,Nx,Ny);

figure(1)
surf(xmat,ymat,TTT'-273)
title('Temperature along the surface')
xlabel('x [m]')
ylabel('y [m]')
zlabel('Temperature [°C]')

figure(2)
contourf(xmat,ymat,TTT'-273,20)
title('Temperature on the surface - projection')
xlabel('x [m]')
ylabel('y [m]')


% Heat power dissipated

qin = trapz(xx,abs(hh*(TT(Ntot-Nx+1:Ntot)-Tair)));
fprintf('The heat power dissipated is %.5f W/m^2', qin)
