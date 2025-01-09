% Steel spheres with a diameter of 12 mm (conductivity = 40 W/m/K, density = 7800 kg/m³,
% and specific heat = 600 J/kg/K) are tempered by rapidly heating them to 1150 K and then
% cooling them slowly to a final temperature of 400 K in an environment with air.
% The air temperature increases over time as Tair = 325 K + 0.0375 K/s × t, where
% t is the time since the start of the cooling process. Numerically calculate the temperature
% evolution of the spheres over time using the implicit Euler scheme, assuming that the
% heat transfer coefficient is 20 W/m²K and that the entire external surface of the spheres
% is in contact with the air. 
% After plotting the time evolution of the sphere temperature, determine how long it takes
% for the sphere temperature to equal that of the air.
% Verify the result obtained using the explicit Euler scheme.

% Laib 8, exercise 5
% Simone Canevarolo
% S269893
% 29/11/2024

clear all
close all
clc

dd = 12e-3; % m
kk = 40; % W/m/K
rovol = 7800; % kg/m^3
cp = 600; % J/kg/K
T0 = 1150; % K
Tend = 400; % K
hh = 20; % W/m^2/K

Tair_funz = @(t) 325+0.0375*t;

% I find the Biot number in order to know if I can solve the problem as a
% zero-dimensional

Bi = hh*(dd/2)/kk;

% Since Biot is Bi<0.1, we can proceed this way.

tmax = 1500; % tempo massimo, arbitrario
tt_air = 1:tmax;

figure(1)
plot(tt_air,Tair_funz(tt_air))
hold on

dt = 1; % s
ii = 1;

TT = T0*ones(tmax,1);
TT(1) = T0;
Tm = T0;

As = 4*pi*dd^2/4; % , area for heat exchange, m^2
Vol = 4/3*pi*(dd/2)^3; % volume, m^3

aa = hh*As/Vol*dt/rovol/cp;
Tair = Tair_funz(TT(1));

while TT(ii)>Tair && ii < tmax

    ii = ii+1;

    Tair = Tair_funz(ii);
    TT(ii) = (Tm+Tair*aa)/(1+aa);
    Tsfera = TT(ii);

    Tm = TT(ii);

end

figure(1)
plot(tt_air,TT)
title('Temperature of the sphere and air vs time')
xlabel('time [s]')
ylabel('Temperature [K]')
legend('Air','Steel')

if TT(end) > Tair
    fprintf('\nThe time required for the sphere temperature to equal that of the air is %d s (BE).\n', ...
        tt_air(end));
    fprintf('\nHowever, the cooling is insufficient as the temperature reached is %.1f K instead of %.1f K.\n', ...
        TT(end), Tair);
else
    fprintf('\nThe time required for the sphere temperature to reach the final temperature is %d s (BE).\n', ...
        tt_air(end));
    fprintf('\nHowever, the cooling is incomplete as the air is still at a temperature of %.1f K.\n', ...
        TT(end));
end
