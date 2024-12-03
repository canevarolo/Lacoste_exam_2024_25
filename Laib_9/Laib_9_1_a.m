% Consider a cylindrical rod of length L = 10 cm and diffusivity alpha = 0.1 cm²/s,
% laterally insulated, where heat propagates only by conduction in the axial
% direction x without any heat generation. At the initial time tin = 0 s,
% the temperature distribution in the rod is
% 
% T(x, 0) = T0 + 50 sin(pi x / L)

% (the domain is assumed from 0 to L), while the ends are maintained at a constant temperature T0 = 300 K.
% Numerically calculate the evolution of the temperature distribution T(x, t)
% along the rod up to the time tend = 100 s.
% 
% Use the finite difference method with the implicit Euler backward scheme,
% with delta t = 0.1 s and delta x = 0.1 cm, and compare the calculated solution
% with the analytical one 
%
% T_an = @(x,t) T0+50*sin(pi*x/ll)*exp(-pi^2/ll^2*alfa*t)
%
% at various time instants (t = 20 s, 40 s, 60 s, 80 s, 100 s).

% Laib 9, exercise 1 A
% Simone Canevarolo
% S269893
% 3/12/2024

clear all
close all
clc

ll = 10e-2; % length, m
alfa = 0.1e-4; % thermal diffusivity, m^2/s
tzero = 0; % time to begin, s

T0 = 300; % temperature of both borders, K
tend = 100; % time to end, s

dt = 0.1; % time step, s
dx = 0.1e-2; % length step, m

T = @(x) T0+50*sin(pi*x/ll); % temperature function

T_an = @(x,t) T0+50*sin(pi*x/ll)*exp(-pi^2/ll^2*alfa*t); % Analytical function of temperature

xx = (0:dx:ll)';
Nx = length(xx);

figure(1)
plot(xx,T(xx),'-*b')
hold on

aa = alfa*dt/dx^2; % coefficient

sup_diag = -aa*ones(Nx,1);
main_diag = (1+2*aa)*ones(Nx,1);
sub_diag = sup_diag;

Band = [[sup_diag(2:end);0], main_diag, [0;sub_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

time = (tzero:dt:tend);

tplot = [20.0 40.0 60.0 80.0 100.0];
tinput = round(tplot/dt)+1;
zz = 1;

Tm = T(xx);
tt = tzero;

for ii = 2:length(time)

    tt = tt+dt;

    bb = Tm;

    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = T0;

    AA(end,end-1) = 0;
    AA(end,end) = 1;
    bb(end) = T0;

    TT = AA\bb;

    if ii == tinput(zz)

        plot(xx,TT,'LineWidth',2)
        hold on

        Tanalytic = T_an(xx,time(ii));
        plot(xx,Tanalytic,'Linewidth',2)
        hold on

        zz = zz+1;

    end

    Tm = TT;

end

legend('0 s','20 s','40 s','60 s','80 s','100 s');
xlabel('Length [m]')
ylabel('Temperature [K]')
title('Temperature at different time')

% Comments:
% we plotted both the analytic function and the function obtained using the
% implicit method of Backward Euler (BE) with precise step of time and space.
% From the plot you can see how the maximum of the temperature is in the middle
% since there is a Dirichlet condition, so fixed temperature, at the edges of the bar.
% There is no visible distinction between the analytic and non-analytic
% function.
