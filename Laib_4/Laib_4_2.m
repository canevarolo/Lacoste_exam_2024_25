% The outer walls of a building consist of a layer of bricks (thermal conductivity kb = 0.55 W/m/K)
% externally coated with a layer of insulation (thermal conductivity ki = 0.04 W/m/K),
% with thicknesses of tb = 25 cm and ti = 10 cm, respectively.
% Calculate and graph the 1D temperature distribution through the wall during winter,
% when the inner wall surface temperature is Tin = 20°C, and the outer surface is exposed
% to convection with air at Text = −5°C (heat transfer coefficient hext = 25 W/m²K).
% Additionally, calculate (and plot) the heat loss through the wall as a function of the
% insulation layer's conductivity, varying it from the value ki to kb.

% Laib 4 , exercise 2
% Simone Canevarolo
% S269893
% 5/11/2024

clear all
close all
clc

llb = 25e-2; % m
lli = 10e-2; % m
lltot = lli+llb; % m
kb = 0.55; % W/m/K
ki = 0.04; % W/m/K
Text = 268; % K
Tin = 293; % K
hh = 25; % W/m^2/K

dx = 1e-3; % m
xx = (0:dx:lltot)';
Nx = length(xx);
 
xxi = (0:dx:lli)';
Nxi = length(xxi);

% xxb = (lli:dx:lltot)';
% Nxb = length(xxb);

% brick
sub_diag = ones(Nx,1);
main_diag = -2*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = zeros(Nx,1);

AA(1,1) = 1+hh*dx/ki;
AA(1,2) = -1;
bb(1) = Text*hh*dx/ki;

AA(Nxi,Nxi-1) = ki;
AA(Nxi,Nxi) = -ki-kb;
AA(Nxi,Nxi+1) = kb;

AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Tin;

TT = AA\bb;

figure(1)
plot(xx,TT-273,'k','LineWidth',2)
hold on
xline(0.1,'r') 
hold off

title('Temperature along the wall')
xlabel('Thickness [m]')
ylabel('Temperature [K]')


%%

kisol = linspace(ki,kb,50);
Nis = length(kisol);
qlosses = zeros(Nis,1);

for ii = 1:Nis

ki = kisol(ii);

dx = 1e-3;
xx = (0:dx:lltot)';
Nx = length(xx);
 
xxis = (0:dx:lli)';
Nxis = length(xxis);

% brick
sub_diag = ones(Nx,1);
main_diag = -2*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = zeros(Nx,1);

AA(1,1) = 1+hh*dx/ki;
AA(1,2) = -1;
bb(1) = Text*hh*dx/ki;

AA(Nxi,Nxi-1) = ki;
AA(Nxi,Nxi) = -ki-kb;
AA(Nxi,Nxi+1) = kb;

AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Tin;

TT = AA\bb;

qlosses(ii) = abs(hh*(Text-TT(1))); 

end

figure(2)
plot(kisol,qlosses,'r-*','LineWidth',2)
title('Heat dispersion vs layer conductivity')
xlabel('Thickness [m]')
ylabel('Heat flux [W/m^2]')