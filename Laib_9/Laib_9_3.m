% A superconducting filament (diameter D = 0.81 mm, length L = 1 m, diffusivity alpha = 11e-5 m^2/s,
% thermal conductivity k = 350 W/m/K), initially at 4.5 K, carries a current I = 50 A.
% The electric field generated by the current while the filament is in a superconducting state is
% E = E0 (I/IC)^n
% where E0 = 1e-5 V/m and n = 1.2.
% IC is the critical current (i.e., the maximum current that can be carried by the conductor
% in the superconducting state), which decreases as the temperature increases: IC = 50 - T^2.
% The ends of the filament are maintained at 4.5 K. Numerically calculate the temperature
% rise up to t = 15000 s using the backward Euler method.
% Plot the time evolution of the volumetric power and the temperature evaluated at the midpoint of the filament.

% Laib 9, exercise 3
% Simone Canevarolo
% S269893
% 28/12/2024

clear all
close all
clc

dd = 0.81e-3; % m
ll = 1; % m
alpha = 11e-5; % m^2/s
kk = 350; % W/m/K
T0 = 4.5; % K

II = 50; % A
nn = 1.2;
E0 = 1e-5; % V/m

Abase = dd^2/4*pi;

% IC = @(T) 50-T^2;
% E = @(IC) E0*(II/IC)^nn;

qv = @(T) (E0*(II/(50-T.^2)).^nn).*II/Abase;

tend = 15000; % s

dx = 1e-3; % m
xx = (0:dx:ll)';
Nx = length(xx);

dt = 1; % s
tt = (0:dt:tend);
Nt = length(tt);

aa = alpha*dt/dx^2;

Tm = T0*ones(Nx,1);
Tmiddle = Tm;
Tmed = T0;

qv_vett = qv(Tmed)*ones(Nt,1);

sub_diag = -aa*ones(Nx,1);
main_diag = (1+2*aa)*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);


for ii = 2:Nt

    bb = Tm+qv(Tmed)*dt*alpha/kk;

    qv_vett(ii) = qv(Tmed);

    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = T0;
    
    AA(end,end-1) = 0;
    AA(end,end) = 1;
    bb(end) = T0;
    
    TT = AA\bb;
    
    Tmed = TT(round(Nx/2)+1);
    
    Tmiddle(ii) = Tmed;
    
    Tm = TT;

end


figure(1)
plot(tt,qv_vett,'b','LineWidth',2)
title('Volumetric heat power vs time')
xlabel('time [s]')
ylabel('Volumetric power [W/m^3]')
grid on

figure(2)
plot(tt,Tmiddle,'r','LineWidth',2)
title('Temperature in the middle point of the bar vs time')
xlabel('time [s]')
ylabel('Temperature [K]')
grid on
