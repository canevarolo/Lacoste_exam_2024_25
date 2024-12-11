% Datas from Laib_9_1_a
% 
% Verify the time convergence of the algorithm (fix delta x at 0.01 m and/or
% 0.005 m and vary delta t across different orders of magnitude)
% by calculating the relative error at time tend

% Laib 9, exercise 1 B
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

dxvett = [0.01 0.005]; % length step, m

dtvett = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20]; % time step, s

T = @(x) T0+50*sin(pi*x/ll); % temperature function

T_an = @(x,t) T0+50*sin(pi*x/ll)*exp(-pi^2/ll^2*alfa*t); % Analytical function of temperature

for yy = 1:length(dxvett)

    dx = dxvett(yy);
    xx = (0:dx:ll)';
    Nx = length(xx);

    Tanalytic = T_an(xx,tend);

for zz = 1:length(dtvett)

    dt = dtvett(zz);

    time = (tzero:dt:tend);
    Nt = length(time);

    tt = tzero;

    Tm = T(xx);

    for ii = 2:Nt
    
        tt = tt+dt;

        aa = alfa*dt/dx^2; % coefficient

        sup_diag = -aa*ones(Nx,1);
        main_diag = (1+2*aa)*ones(Nx,1);
        sub_diag = sup_diag;
        
        Band = [[sup_diag(2:end);0], main_diag, [0;sub_diag(1:end-1)]];
        
        AA = spdiags(Band,-1:1,Nx,Nx);
    
        bb = Tm;
    
        AA(1,1) = 1;
        AA(1,2) = 0;
        bb(1) = T0;
    
        AA(end,end-1) = 0;
        AA(end,end) = 1;
        bb(end) = T0;
    
        TT = AA\bb;
    
        Tm = TT;
    
    end

    err(zz) = norm(TT-Tanalytic)/norm(Tanalytic-T0);

end

figure(1)
loglog(dtvett,err,'LineWidth',2)
hold on

end

legend('0.01 m', '0.005 m','Location','best');
xlabel('Time step [s]')
ylabel('Relative error [-]')
title('Time convergence at different space steps')
