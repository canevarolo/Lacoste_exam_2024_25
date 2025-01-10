% A rod with a diameter of D = 10 mm and a length of L = 100 mm has one end maintained at Tb = 100 °C.
% The surface of the rod is cooled by free convection with ambient air at Tinf = 25°C and a heat transfer
% coefficient that depends on the temperature difference between the surface and the ambient air according
% to the correlation:  
% 
% hfc = 2.89 * [0.6 + 0.624 * (T - Tinf)^(1/6)]^2  
% 
% where the units are [W/m²/K] for the heat transfer coefficient and [K] for the temperatures.  
% The surface of the rod has an emissivity of ε = 0.2 and exchanges heat radiatively with the surroundings
% at a fixed temperature of Tsurf = 25°C.  
% Using the finite difference method in 2D, calculate:  
% a. The temperature distribution.  
% b. The average temperature at the tip.

% Laib 10, exercise 7
% Simone Canevarolo
% S269893
% 10/01/2025

clear all
close all
clc

dd = 10e-3; % m
ll = 100e-3; % m

Tb = 373; % K
Tair = 298; % K

hh = @(T) 2.89 * [0.6 + 0.624 * (T - Tair)^(1/6)]^2 ;

emiss = 0.2;

Tsur = 298; % K
cond = 14; % W/m/K

dx = 1e-3; % m
dy = dx;

xx = (0:dx:dd)';
yy = (0:dy:ll)';

Nx = length(xx);
Ny = length(yy);

Ntot = Nx*Ny;

[xmat,ymat] = meshgrid(xx,yy);

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2/dx^2-2/dy^2;
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

bb = zeros(Ntot,1);
Tm = Tb*ones(Ntot,1);

tend = 1000;

dt = 1; % s
tt = (0:dt:tend);
Nt = length(tt);

for ii = 2:Nt
    
    % West
    
    kk = 1:Nx;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Tb;
    
    % South
    
    ii = Nx;
    for jj = 2:Ny-1
    
        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = cond/2/dy^2;
        AA(kk,kk-1) = cond/dx^2;
        AA(kk,kk) = -cond/dx^2-cond/dy^2-hh(Tm(kk))/dx;
        AA(kk,kk+Nx) = cond/2/dy^2;
        bb(kk) = -Tair*hh(Tm(kk))/dx;
    
    end
    
    
    % East

    %%% CONDIZIONE IRRAGGIAMENTO DA INSERIRE %%%

    
    
    % North
    
    ii = 1;
    for jj = 2:Ny-1
    
        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = cond/2/dy^2;
        AA(kk,kk) = -cond/dx^2-cond/dy^2-hh(Tm(kk))/dx;
        AA(kk,kk+1) = cond/dx^2;
        AA(kk,kk+Nx) = cond/2/dy^2;
        bb(kk) = -Tair*hh(Tm(kk))/dx;
    
    end
    
    
    TT = AA\bb;
    
    Tm = TT;

end

TTT = reshape(TT,Nx,Ny);

figure(1)
surf(xmat,ymat,TTT'-273)
