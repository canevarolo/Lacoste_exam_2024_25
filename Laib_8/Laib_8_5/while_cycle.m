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
% 09/01/2025

dd = 12e-3; % m
kk = 40; % W/m/K
rovol = 7800; % kg/m^3
cp = 600; % J/kg/K
hh = 20; % W/m^2/K

As = 4*pi/4*dd^2; % m^2
VV = 4/3*pi*(dd/2)^3; % m^3

T0 = 1150;
Tend = 400;
Taria = @(tt) 325+0.0375*tt;

dt = 1;
tt = 0:dt:1500;
mm = length(tt);

Tm = T0;%Tm
TT = Tm*ones(mm,1); %Tm+1

aa = hh*As*dt/VV/rovol/cp;

toll = 1e-7;
err = 10*toll;
ii = 1;
itermax = 1e4;

while err>toll & ii<itermax

    ii = ii+1;
    TT(ii) = (Tm+aa*Taria(ii))/(1+aa);
    tt(ii) = tt(ii-1)+dt;
    err(ii) = norm(TT(ii)-Tm)/norm(TT(ii));
    Tm = TT(ii);

end

figure(1)
plot(tt,TT,'LineWidth',2)
hold on
plot(tt,Taria(tt),'LineWidth',2)

if ii < itermax
    fprintf('The time required is %d seconds', ii)
else
    fprintf('The time required is more than %d seconds', itermax)
end
