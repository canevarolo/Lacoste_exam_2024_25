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

%%

ll = 0.2; % m
base = 0.1; % m
alt = 0.1; % m
alt1 = 0.05; %  m

roacc = 7800; % kg/m^3
cond = 1; % W/m/K
cp = 450; % J/kg/K

Tout = 5+273; % K
Tin = 20+273; % K

%%

dx = 1e-2; % m
dy = dx;

xx = (0:dx:alt)';
xx1 = (0:dx:alt1)';

yy = (0:dy:ll)';
yy1 = (0:dy:base)';

Nx = length(xx);
Nx1 = length(xx1);

Ny = length(yy);
Ny1 = length(yy1);

Ntot = Nx*Ny1+Nx1*(Ny-Ny1);
Ntot1 = Nx*Ny1;
Ntot2 = Nx1*(Ny-Ny1);

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);


bb = zeros(Ntot,1);

% Left rectangle

for ii = 2:Nx-1
    for jj = 2:Ny1-1

        kk = Nx*(jj-1)+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

% Left rectangles

for ii = 2:(Nx1-1)
    for jj = Ny1:Ny-1

        if jj == Ny1   

            kk = Ntot1-Nx+ii;
            AA(kk,kk-Nx) = 1/dy^2;
            AA(kk,kk-1) = 1/dx^2;
            AA(kk,kk) = -2*(1/dx^2+1/dy^2);
            AA(kk,kk+1) = 1/dx^2;
            AA(kk,kk+Nx) = 1/dy^2;

        elseif  jj == Ny1+1
               
            kk=Ntot1+ii;
            AA(kk,kk-Nx)=1/dy^2;
            AA(kk,kk-1)=1/dx^2;
            AA(kk,kk)=-2*(1/dx^2+1/dy^2);
            AA(kk,kk+1)=1/dx^2;
            AA(kk,kk+Nx1)=1/dy^2;
                   
        else

            kk = Ntot1+(jj-Ny1-1)*Nx1+ii;
            AA(kk,kk-Nx1) = 1/dy^2;  
            AA(kk,kk-1) = 1/dx^2;
            AA(kk,kk) = -2*(1/dx^2+1/dy^2);
            AA(kk,kk+1) = 1/dx^2;
            AA(kk,kk+Nx1) = 1/dy^2;

        end
    end
end


% West (left square) - Dirichlet condition
jj = 1;
for ii = 2:Nx

    kk = ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Tout;

end

% North (left square), Dirichlet condition

ii = 1;
for jj = 1:Ny1

    kk = (jj-1)*Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Tout;

end

% South-East (left square) - Dirichlet condition

jj = Ny1;
for ii = Nx1:Nx

    kk = (jj-1)*Nx+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Tin;

end

% South (left square) - Neumann condition

ii = Nx;
for jj = 2:Ny1-1

    kk = (jj-1)*Nx+ii;
    AA(kk,kk-Nx) = 1/dy^2;
    AA(kk,kk-1) = 2*(1/dx^2);
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+Nx) = 1/dy^2;

end

% North (right rectangle) - Dirichlet condition

ii=1;
for jj = Ny1+1:Ny

    kk = Ntot1+(jj-Ny1-1)*Nx1+ii;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Tout;

end

% East (right rectangle) - Neumann adiabatico

jj = Ny;
for ii = 2:Nx1-1

    kk = Ntot-Nx1+ii;
    AA(kk,kk-Nx1) = 2*(1/dy^2);
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+1) = 1/dx^2;

end

% South (right rectangle) - Dirichlet condition

ii = Nx1;
for jj = Ny1+1:Ny

     kk = Ntot1+(jj-Ny1-1)*Nx1+ii;
     AA(kk,kk) = eye(length(kk));
     bb(kk) = Tin;

end

%%

TT = AA\bb;

TTT1 = reshape(TT(1:Ntot1),Nx,Ny1);
TTT2 = reshape(TT(Ntot1+1:end),Nx1,Ny-Ny1);
TTT3 = NaN*ones(Nx-Nx1,Ny-Ny1);

TTT = [TTT1,[TTT2;TTT3]];

figure (1)
surf(xmat,ymat,TTT'-273)

xlabel('x [m]','fontsize',16)
ylabel('y [m]','fontsize',16)
zlabel('Temperature [°C]','fontsize',16)
title('Temperature along the surface')
