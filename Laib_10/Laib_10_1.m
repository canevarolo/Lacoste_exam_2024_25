% An industrial furnace is supported by a brick column with a cross-section of 1 m × 1 m.
% During steady-state operation, the furnace setup ensures that three faces of the column are maintained at 500 K,
% while the remaining face is exposed to an air stream at 300 K, with a heat transfer coefficient of h = 10 W/m²K.
% Determine the temperature map across the column's cross-section and the heat power dissipated into the air per unit length of the column.
% Density = 1920 kg/m^3; Thermal conductivity = 0.72 W/m/K; Specific heat = 835 J/kg/K

% Laib 10, exercise 1
% Simone Canevarolo
% 31/12/2024
% S269893

clear all
close all
clc

ll = 1; % m
Tlato = 500; % K
Tair = 300; % K
hh = 10; % W/m^2/K

rovol = 1920; % Kg/m^3
cond = 0.72; % W/m/K
cp = 835; % J/kg/K

dx = 5e-2; % m
dy = dx;

xx = (0:dx:ll)';
yy = xx;

Nx = length(xx);
Ny = Nx;
Ntot = Nx*Ny;

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[], Ntot, Ntot, 5*Ntot);

bb = zeros(Ntot,1);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = Nx*(jj-1)+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end


% West, mantained temperature (Dirichlet condition)
kk = 1:Nx;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tlato;

% South, mantained temperature (Dirichlet condition)
kk = Nx:Nx:Ntot;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tlato;

% East, mantained temperature (Dirichlet condition)
kk = Nx*(Ny-1)+1:Ntot;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Tlato;

% North, air convection (Robin condition)
ii = 1;
for jj = 2:Ny-1

    kk = Nx*(jj-1)+ii;
    AA(kk,kk-Nx) = cond/2/dy^2;
    AA(kk,kk) = -(cond*(1/dx^2+1/dy^2)+hh/dx);
    AA(kk,kk+1) = cond/dx^2;
    AA(kk,kk+Nx) = cond/2/dy^2;
    bb(kk) = -Tair*hh/dx;

end


TT = AA\bb;
TTT = reshape(TT,Nx,Ny);

figure(1)
surf(xmat',ymat',TTT)
colorbar
xlabel('x [m]')
ylabel('y [m]')
zlabel('Temperature [K]')
title('Temperature on the section')



deltax = [dx/2, dx*ones(1,Nx-2), dx/2];
deltat = (TTT(:,Ny-1)-TTT(:,Ny));
qexchange = kk*deltat/dy;
qtot = sum(qexchange.*deltax,"all");

fprintf("Heat exchange is %.5f W/m^2", qtot)
