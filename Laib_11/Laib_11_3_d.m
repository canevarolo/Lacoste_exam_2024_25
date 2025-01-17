% d. Calculate the temperature map of the bar in the steady state reached after the heater is turned on,
% identifying the thermal level and the location of the hottest point. Solve this both as a transient problem
% and as a steady-state problem with the heater on.

% Laib 11, exercise 1
% Simone Canevarolo
% S269893
% 17/01/2025

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

dx = 5e-3; % m
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

dt = 1; % s
tend = 1000; % s
tt = (0:dt:tend);

toll = 1e-4;
err = 1+toll;

Tm = Th2o*ones(Ntot,1);

while err>toll

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


    % North East - heat generation + Robin

    ii = 1;
    for jj = Ny1:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = cond/2/dy^2;
        AA(kk,kk) = -cond/dx^2-cond/dy^2-hh/dx;
        AA(kk,kk+1) = cond/dx^2;
        AA(kk,kk+Nx) = cond/2/dy^2;
        bb(kk) = -Taria*hh/dx-qvol/dx;

    end


    % Northeast - Robin + Robin + heat generation
    
    kk = Ntot-Nx+1;
    AA(kk,kk-Nx) = cond/2/dy^2;
    AA(kk,kk) = -cond/2/dx^2-cond/2/dy^2-hh*(1/dx+1/dy);
    AA(kk,kk+1) = cond/2/dx^2;
    bb(kk) = -Taria*hh*(1/dx+1/dy)-qvol/dy/2;
    
    
    TT = AA\bb;
    TTT = reshape(TT,Nx,Ny);

    err = norm(TT-Tm)/norm(TT-Taria);


    figure(1)
    surf(xmat,ymat,TTT'-273)

    Tm = TT;

end


figure(2)
contourf(xmat,ymat,TTT'-273,20)
