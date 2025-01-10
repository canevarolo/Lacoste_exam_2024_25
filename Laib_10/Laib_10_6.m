% A straight fin with a constant cross-section is made of a material with thermal conductivity k = 50 W/m/K,
% thickness w = 6 mm, and length L = 48 mm. It is very long in the direction perpendicular to the figure.
% The convective heat transfer coefficient is h = 500 W/m²/K, and the ambient air temperature is Ta = 30°C.
% The base of the fin is maintained at a constant temperature of Tb = 100°C, while the tip is perfectly insulated.

% a. Using the finite difference method, estimate the temperature distribution along the fin.
% Is the assumption of one-dimensional heat transfer reasonable for this fin?

% b. Estimate the heat dissipated by the fin per unit length normal to the figure [W/m].

% c. Using the mesh from part (a), calculate and plot the temperature distribution for h = 10, 100, 500, and 1000 W/m²/K.

% Lab 10, exercise 5
% Simone Canevarolo
% S269893
% 09/01/2025

clear all
close all
clc

spess = 6e-3; % m
lungh = 48e-3; % m
cond = 50; % W/m/K
Taria = 30 + 273; % K
Tbase = 100 + 273; % K

hvett = [10, 100, 500, 1000]; % W/m²/K

dx = 1e-3; % m
dy = dx;

xx = (0:dx:spess)';
yy = (0:dy:lungh)';

Nx = length(xx);
Ny = length(yy);

Ntot = Nx * Ny;

[xmat, ymat] = meshgrid(xx, yy);

AA = sparse([], [], [], Ntot, Ntot, 5 * Ntot);

for zz = 1:length(hvett)

    for ii = 2:Nx-1
        for jj = 2:Ny-1
    
            kk = (jj-1) * Nx + ii;
            AA(kk, kk-Nx) = 1 / dy^2;
            AA(kk, kk-1) = 1 / dx^2;
            AA(kk, kk) = -2 * (1 / dx^2 + 1 / dy^2);
            AA(kk, kk+1) = 1 / dx^2;
            AA(kk, kk+Nx) = 1 / dy^2;
    
        end
    end
    
    bb = zeros(Ntot, 1);
    
    % West boundary
    kk = 1:Nx;
    AA(kk, :) = 0;
    AA(kk, kk) = eye(length(kk));
    bb(kk) = Tbase;
    
    % South boundary
    ii = Nx;
    for jj = 2:Ny-1

        kk = (jj-1) * Nx + ii;
        AA(kk, kk-Nx) = 1 / (2 * dy^2);
        AA(kk, kk-1) = 1 / dx^2;
        AA(kk, kk) = -1 / dx^2 - 1 / dy^2 - hvett(zz) / cond / dx;
        AA(kk, kk+Nx) = 1 / (2 * dy^2);
        bb(kk) = -Taria * hvett(zz) / cond / dx;
        
    end
    
    % North boundary
    jj = Ny;
    for ii = 2:Nx-1

        kk = (jj-1) * Nx + ii;
        AA(kk, kk-Nx) = 1 / dy^2;
        AA(kk, kk-1) = 1 / (2 * dx^2);
        AA(kk, kk) = -1 / dx^2 - 1 / dy^2;
        AA(kk, kk+1) = 1 / (2 * dx^2);

    end
    
    % East boundary
    ii = 1;
    for jj = 2:Ny-1

        kk = (jj-1) * Nx + ii;
        AA(kk, kk-Nx) = 1 / (2 * dy^2);
        AA(kk, kk) = -1 / dx^2 - 1 / dy^2 - hvett(zz) / cond / dx;
        AA(kk, kk+1) = 1 / dx^2;
        AA(kk, kk+Nx) = 1 / (2 * dy^2);
        bb(kk) = -Taria * hvett(zz) / cond / dx;

    end
    
    % North-East corner
    kk = Ntot - Nx + 1;
    AA(kk, kk-Nx) = cond / dy^2;
    AA(kk, kk) = -cond / dx^2 - cond / dy^2 - hvett(zz) / dx^2;
    AA(kk, kk+1) = cond / dx^2;
    bb(Ntot - Nx + 1) = -Taria * hvett(zz) / dx^2;
    
    % South-East corner
    kk = Ntot;
    AA(kk, kk-Nx) = cond / dx^2;
    AA(kk, kk-1) = cond / dy^2;
    AA(kk, kk) = -cond / dx^2 - cond / dy^2 - hvett(zz) / dx^2;
    bb(kk) = -Taria * hvett(zz) / dx^2;
    
    TT = AA \ bb;
    
    TTT = reshape(TT, Nx, Ny);

    if hvett(zz) == 500

        figure(1)
        surf(xmat, ymat, TTT' - 273)
        title('Temperature distribution on the plate at h = 500 W/m²/K')
        xlabel('X [m]')
        ylabel('Y [m]')
        zlabel('Temperature [K]')
            
        figure(2)
        contourf(xmat, ymat, TTT' - 273, 50)
        title('2D temperature plot of the plate at h = 500 W/m²/K')
        xlabel('X [m]')
        ylabel('Y [m]')

    end

    figure(3)
    plot(yy, TT(1:Nx:end) - 273)
    hold on

    qout = abs(trapz(yy, hvett(zz) * (TT(1:Nx:end) - Taria)));
    fprintf('Heat dissipated with heat transfer coefficient of %d W/m^2/K is %.5f W/m²\n', hvett(zz), qout)

end

legend('10 W/m²/K', '100 W/m²/K', '500 W/m²/K', '1000 W/m²/K')
title('Temperature profile on the piece for different h values')
xlabel('Length [m]')
ylabel('Temperature [K]')


```matlab
% 1D assumption is reasonable; the projection of the T graph along the length
% appears almost one-dimensional.
```
