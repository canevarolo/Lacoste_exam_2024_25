% A very long bar made of a material with a density of 170 kg/m³, specific heat of 210 J/kgK,
% and thermal conductivity of 2.5 W/mK has a rectangular cross-section and floats on the free surface
% of a large water basin.
% Initially, the bar is in thermal equilibrium with the surrounding environment (air at a temperature of 30°C),
% exchanging heat on all sides above the water surface with a heat transfer coefficient of 10 W/m²K,
% and is perfectly thermally coupled with the water at 5°C. 
% At a certain moment, a heater with negligible thickness is activated on the upper edge of the bar,
% in the position shown in Figure 2. The heater supplies the bar with a heat flux of 1 kW/m²,
% initiating a transient thermal process that alters the temperature distribution within the bar. 
% The thermal behavior of the bar is to be studied, specifically including the following on the generic 2D section:
% 
% a. Calculate the temperature map of the bar in the initial steady-state condition, identifying the thermal level
% and the location of the hottest point.
% 
% b. Perform an engineering study of the accuracy of the thermal level of the hottest point in the bar in the
% initial steady state, varying the spatial discretization parameters.
% 
% c. Analyze the thermal balance of the bar in the initial steady state to verify that the power entering the bar
% through the sides exposed to convection with the air is balanced by the power dissipated through the sides of the bar
% in contact with water.


% Laib 11, exercise 1
% Simone Canevarolo
% S269893
% 16/01/2025

clear all
close all
clc

base = 25e-2; % m
alt = 18e-2; % m
alt1 = 15e-2; % m
base1 = 20e-2; % m

rovol = 170; % kg/m^3
cp = 210; % K/kg/K
cond = 2.5; % W/m/K
Taria = 30+273; % K
hh = 10; % W/m^2/K
Th2o = 5+273; % K

qvol = 1e3; % W/m^3

dxvett = [1e-2 5e-3 2e-3 1e-3];
Tmaxvett = zeros(length(dxvett),1);

for zz = 1:length(dxvett)

    dx = dxvett(zz);
    dy = dx;

    xx = (0:dx:alt)';
    yy = (0:dy:base)';
    
    xx1 = (0:dx:alt1)';
    yy1 = (0:dy:base1)';
    
    Nx = length(xx);
    Ny = length(yy);
    
    Nx1 = length(xx1);
    Ny1 = length(yy1);
    
    Ntot = Nx*Ny;
    
    AA = sparse([],[],[],Ntot,Ntot,5*Ntot);
    [xmat,ymat] = meshgrid(xx,yy);


    for ii = 2:Nx-1
        for jj = 2:Ny-1
    
            kk = (jj-1)*Nx+ii;
            AA(kk,kk-Nx) = 1/2/dy^2;
            AA(kk,kk-1) = 1/2/dx^2;
            AA(kk,kk) = -1/dx^2-1/dy^2;
            AA(kk,kk+1) = 1/2/dx^2;
            AA(kk,kk+Nx) = 1/2/dy^2;
    
        end
    end
    
    bb = zeros(Ntot,1);
    
    
    % West - Robin
    
    jj = 1;
    for ii = 2:Nx1-1
    
        kk = (jj-1)*Nx+ii;
        AA(kk,kk-1) = cond/2/dx^2;
        AA(kk,kk) = -cond/dx^2-cond/dy^2-hh/dy;
        AA(kk,kk+1) = cond/2/dx^2;
        AA(kk,kk+Nx) = cond/dy^2;
        bb(kk) = -Taria*hh/dy;
    
    end
    
    
    % North - Robin
    
    ii = 1;
    for jj = 2:Ny-1
    
        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = cond/2/dy^2;
        AA(kk,kk) = -cond/dx^2-cond/dy^2-hh/dx;
        AA(kk,kk+1) = cond/dx^2;
        AA(kk,kk+Nx) = cond/2/dy^2;
        bb(kk) = -Taria*hh/dx;
    
    end
    
    
    % East - Robin
    
    jj = Ny;
    for ii = 2:Nx1-1
    
        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = cond/dy^2;
        AA(kk,kk-1) = cond/2/dx^2;
        AA(kk,kk) = -cond/dx^2-cond/dy^2-hh/dy;
        AA(kk,kk+1) = cond/2/dx^2;
        bb(kk) = -Taria*hh/dy;
    
    end
    
    
    % Northwest - Robin + Robin
    
    kk = 1;
    AA(kk,kk) = -cond/2/dx^2-cond/2/dy^2-hh*(1/dx+1/dy);
    AA(kk,kk+1) = cond/2/dx^2;
    AA(kk,kk+Nx) = cond/2/dy^2;
    bb(kk) = -Taria*hh*(1/dx+1/dy);
    
    
    % Northeast - Robin + Robin
    
    kk = Ntot-Nx+1;
    AA(kk,kk-Nx) = cond/2/dy^2;
    AA(kk,kk) = -cond/2/dx^2-cond/2/dy^2-hh*(1/dx+1/dy);
    AA(kk,kk+1) = cond/2/dx^2;
    bb(kk) = -Taria*hh*(1/dx+1/dy);
    
    
    % West (water) - Dirichlet
    
    kk = Nx1:Nx;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;
    
    
    % South (water) - Dirichlet
    
    kk = Nx:Nx:Ntot;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;
    
    
    % East (water) - Dirichlet
    
    kk = Ntot-Nx+Nx1:Ntot;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    bb(kk) = Th2o;
    
    
    TT = AA\bb;
    
    Tmaxvett(zz) = max(TT);

end

TTT = reshape(TT,Nx,Ny);

% Point A
figure(1)
surf(xmat,ymat,TTT'-273)

figure(2)
contourf(xmat,ymat,TTT'-273,20)

% Point B
figure(3)
loglog(dxvett,Tmaxvett-273)
