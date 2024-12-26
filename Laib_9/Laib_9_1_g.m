% Datas from Laib_9_1_a

% Evaluate the time evolution of the exchanged power and the total energy transferred after 100 seconds.
% Assume the bar is made of steel with a thermal conductivity of 10 W/m/K.

% Laib 9, exercise 1 E
% Simone Canevarolo
% 23/12/2024
% S269893

clear all
close all
clc

ll = 10e-2; % length, m
alfa = 0.1e-4; % thermal diffusivity, m^2/s
tzero = 0; % time to begin, s
kk = 10; % W/m/K

T0 = 300; % temperature of both borders, K
tend = 100; % time to end, s

dt = 0.1; % time step, s
dx = 0.1e-2; % length step, m

T = @(x) T0+50*sin(pi*x/ll); % temperature function

xx = (0:dx:ll)';
Nx = length(xx);

aa = alfa*dt/dx^2; % coefficient

sup_diag = -aa*ones(Nx,1);
main_diag = (1+2*aa)*ones(Nx,1);
sub_diag = sup_diag;

Band = [[sup_diag(2:end);0], main_diag, [0;sub_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

time = (tzero:dt:tend);

Eout = 0; % 

Tm = T(xx);
tt = tzero;

Qoutvett = zeros(length(time),1);

for ii = 1:length(time)

    tt = tt+dt;

    bb = Tm;

    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = T0;

    AA(end,end-1) = 0;
    AA(end,end) = 1;
    bb(end) = T0;

    TT = AA\bb;


    % There is no convection since the bar is laterally insulated.
    % There is no heat generation inside the bar, the power is produced by
    % conduction

    Qoutvett(ii) = abs(kk*(TT(2)-TT(1))/dx)+abs(kk*(TT(end-1)-TT(end))/dx);

    Eout = Eout + Qoutvett(ii)*dt;

    Tm = TT;

end

figure(1)
plot(time,Qoutvett,'r','LineWidth',2)
title('Power exchange vs time')
xlabel('time [s]')
ylabel('Power released [W/m^2]')
legend('Power','Location','best')

% print the final energy output value
Eout