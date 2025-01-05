% Consider the cross-section of a steel beam, shown in Figure 1.

% a. Compute and represent the temperature map under steady-state conditions, 
% with the boundary conditions depicted in the figure and Tin = 20 °C.

% b. Calculate the heat power dissipated per unit length of the beam through the edges 
% that are at a lower temperature.

% For steel, assume a density of 7800 kg/m³, a thermal conductivity of 1 W/mK, 
% and a specific heat capacity of 450 J/kgK.


% Laib 10, exercise 2
% Simone Canevarolo
% S269893
% 04/01/2025

clear all
close all
clc

alt = 0.1; % m
lungh = 0.2; % m
base = 0.1; % m
alt1 = 0.05; % m

rovol = 7800; % kg/m^3
kk = 1; % W/m/K
cp = 450; % J/kg/K

Tin = 20+273; % K
Tout = 5+273; % K

dx = 1e-4; % m
dy = dx;

xx = (0:dx:alt)';
yy = (0:dy:lungh)';
xx1 = (0:dx:alt1)';
yy1 = (0:dy:base)';

Nx = length(xx);
Ny = length(yy);
Nalt1 = length(xx1);
Nybase = length(yy1);

Nleft = Nx*Nybase;
Nrightup = Nleft/2;
Nrightdown = Nrightup;

Ntot = Nx*Ny;

% Left rectangle

for ii = 2:Nybase-1
    for jj = 2:Nx-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

% Right up rectangle

for ii = 2:Nalt1-1
    for jj = Nbase+2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

% Right down rectangle

for ii = Nalt1+2:Nx-1
    for jj = Nbase+2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end


bb = zeros(Ntot,1);


