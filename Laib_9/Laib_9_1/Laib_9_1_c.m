% Datas from Laib_9_1_a
% 
% Verify the spatial convergence of the algorithm (fixing ∆t at 0.1 s
% and/or 0.05 s and varying ∆x across several orders of magnitude)
% by calculating the relative error at time tend.

% Laib 9, exercise 1 C
% Simone Canevarolo
% S269893
% 11/12/2024

clear all
close all
clc

ll = 10e-2; % length, m
alfa = 0.1e-4; % thermal diffusivity, m^2/s
tzero = 0; % time to begin, s

T0 = 300; % temperature of both borders, K
tend = 100; % time to end, s

dtvett = [0.05 0.1]; % time step, m

dxvett = [5e-6 1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2]; % length step, m

Tfunz = @(x) T0+50*sin(pi*x/ll); % temperature function

T_an = @(x,t) T0+50*sin(pi*x/ll)*exp(-pi^2/ll^2*alfa*t); % Analytical function of temperature

for yy = 1:length(dtvett)

    dt = dtvett(yy);
    time = (tzero:dt:tend);
    Nt = length(time);

    for zz = 1:length(dxvett)
    
        dx = dxvett(zz);
        xx = (0:dx:ll)';
        Nx = length(xx);

        Tanalytic = T_an(xx,tend);

        tt = tzero;
    
        Tm = Tfunz(xx);
    
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
    loglog(dxvett,err,'LineWidth',2)
    hold on

end

legend('0.05 s','0.1 s','Location','best');
xlabel('Space step [m]')
ylabel('Relative error [-]')
title('Spacial convergence at different time steps')