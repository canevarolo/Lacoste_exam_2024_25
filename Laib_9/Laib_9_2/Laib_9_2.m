% An electric current of 50 A flows through an electrical conductor
% (laterally insulated, with a length L = 1 m, resistivity r = 1.75 × 10⁻⁸ Ω·m,
% diffusivity a = 11 × 10⁻⁵ m²/s, thermal conductivity k = 350 W/(m·K)),
% with an internal diameter D = 4 mm. Calculate the temperature evolution
% at the midpoint of the conductor, initially at a temperature T(x, 0) = T₀ = 300 K.
% The right end of the conductor is maintained at a constant temperature T₀,
% while the temperature at the left end is measured by a sensor; the measurements
% are available in the file Tleft.dat.
% Solve the problem numerically using finite differences and an implicit time-stepping scheme.

% Laib 9, exercise 2
% Simone Canevarolo
% S269893
% 26/12/2024

clear all
close all
clc

II = 50; % electric current, A
ll = 1; % m
roel = 1.75e-8; % Ohm*m
alpha = 11e-5; %m^-2/s
kk = 360; % W/m/K
dd = 4e-3; % m

T0 = 300; % K

dx = 1e-3; % m
xx = (0:dx:ll)';
Nx = length(xx);

xmed = round(Nx/2)+1;

Tmidvett = T0*ones(Nx,1);
Tm = T0*ones(Nx,1);

qv = II^2*roel*ll/(dd^2/4*pi); % W/m^3
cp = kk/roel/alpha; % J/Kg/K

edge = load('Tleft.dat');
tt = edge(:,1);
Tsxvett = edge(:,2);

dt = 10; % s
Nt = length(tt);

aa = alpha*dt/dx^2;

sub_diag = -aa*ones(Nx,1);
main_diag = (1+2*aa)*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);
    
for ii = 2:Nt

    bb = Tm+qv*dt/roel/cp*ones(Nx,1);

    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = Tsxvett(ii);

    AA(end,end-1) = 0;
    AA(end,end) = 1;
    bb(end) = T0;
       
    TT = AA\bb;

    Tmidvett(ii) = TT(xmed);
    
    Tm = TT;

end

figure(1)
plot(tt,Tmidvett-T0)


% DA TERMINARE!