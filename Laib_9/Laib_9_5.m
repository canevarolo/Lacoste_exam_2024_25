% Steel spheres with a diameter of 12 mm (conductivity = 40 W/m/K, density = 7800 kg/m³, and specific heat = 600 J/kg/K)
% are tempered by rapidly heating them to 1150 K and then cooling them slowly to a final temperature of 400 K
% in an environment with air, whose temperature T_air increases over time: T_air = 325 K + 0.0375 K/s × t,
% where t is the time elapsed since the beginning of the cooling process.
% Calculate numerically, using the implicit Euler scheme, the time evolution of the temperature of the spheres,
% assuming a heat transfer coefficient of 20 W/m²K and that the entire external surface of the spheres is in
% contact with the air.
% Perform a one-dimensional analysis to show the time evolution of the average temperature of the spheres and
% determine after how long it becomes equal to that of the air. Verify the obtained result using the forward Euler scheme.
% Finally, compare the obtained results with the zero-dimensional analysis conducted in Exercise 4 of Lab #7.

% Laib 9, exercise 5
% Simone Canevarolo
% S269893
% 30/12/2024

clear all
close all
clc

dd_tot = 12e-3; % m
dd = dd_tot/2;

kk = 40; % W/m/K
rovol = 7800; % kg/m^3
cp = 600; % J/kg/K

T0 = 1150; % K
Tend = 400; % K
hh = 20; % W/m^2/K

% Temperature of the air over time (t)
T_air = @(t) 325+0.0375*t;

dx = 1e-4; % m
xx = (0:dx:dd)';
Nx = length(xx);

dt = 1; % s
tt = (0:dt:5000);
Nt = length(tt);

Tm = T0*ones(Nx,1);
Tmedia = T0*ones(Nt,1);

aa = kk*dt/rovol/cp/dx^2;

% Solve the problem using BE method (implicit)

for ii = 2:Nt

    sub_diag = -aa*ones(Nx,1);
    main_diag = (1+2*aa)*ones(Nx,1);
    sup_diag = sub_diag;
    
    Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
    
    AA = spdiags(Band,-1:1,Nx,Nx);
    
    bb = Tm;

    Taria = T_air(ii);
    
    AA(1,1) = kk/dx+hh;
    AA(1,2) = -kk/dx;
    bb(1) = hh*Taria;
    
    AA(end,end-1) = -1;
    AA(end,end) = 1;
    bb(end) = 0;
    
    TT = AA\bb;
    
    Tmedia(ii) = mean(TT);
    
    Tm = TT;

end

figure(1)
plot(tt,Tmedia,'LineWidth',2)
hold on
plot(tt,T_air(tt), 'LineWidth',2)
title('Temperature of the steel balls and air vs time')
xlabel('time [s]')
ylabel('Temperature [K]')
legend('Steel balls', ' Air', 'Location','best')



%% TERMINARE CON FORWARD EULER E L'ANALISI 0-DIMENSIONALE!
