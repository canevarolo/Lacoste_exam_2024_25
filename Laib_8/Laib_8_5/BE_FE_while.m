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
% 09/01/2024


clear all
close all
clc

dd = 12e-3; % m
kk = 40; % W/m/K
rovol = 7800; % kg/m^3
cp = 600; % J/kg/K
hh = 20; % W/m^2/K

T0 = 1150; % K
Tend = 400; % K

Taria = @(t) 325+0.0375*t;

dt = 1; % s
tt = (0:dt:1500);
Nt = length(tt);

aa = hh*3/(dd/2)*dt/rovol/cp;

Tm = T0;
TBE = T0*ones(Nt,1);

toll = 1e-7;
err_be = 10*toll;
ii = 1;
itermax = 1e4;

while err_be>toll & ii<itermax

    ii = ii+1;
    TBE(ii) = (Tm+aa*Taria(ii))/(1+aa);
    tt(ii) = tt(ii-1)+dt;
    err_be(ii) = norm(TBE(ii)-Tm)/norm(TBE(ii));
    Tm = TBE(ii);

end

figure(1)
plot(tt,TBE,'LineWidth',2)
hold on
plot(tt,Taria(tt),'LineWidth',2)
hold on


Tm = T0;
TFE = T0*ones(Nt,1);

tt_fe = (0:dt:1500); % s

toll = 1e-7;
err_fe = 10*toll;
ii = 1;
itermax = 1e4;

while err_fe>toll & ii<itermax

    ii = ii+1;
    TFE(ii) = aa*Taria(ii)-Tm*(aa-1);
    tt_fe(ii) = tt_fe(ii-1)+dt;
    err_fe(ii) = norm(TFE(ii)-Tm)/norm(TFE(ii));
    Tm = TFE(ii);

end

plot(tt_fe,TFE,'g--','LineWidth',2)


