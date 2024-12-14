% Datas from 9_1_a

% Evaluate the temperature evolution over time using the Crank-Nicolson
% method and verify its temporal convergence by fixing Δx = 0.001 m.

% Laib 9, exercise 1 E
% Simone Canevarolo
% 14/12/2024
% S269893

clear all
close all
clc

ll = 10e-2; % length, m
alpha = 1e-5; % thermal diffusivity, m^2/s

T0 = 300; % temperature of both borders, K
tend = 100; % time to end, s

T = @(x) T0+50*sin(pi*x/ll); % temperature function

T_an = @(x,t) T0+50*sin(pi*x/ll)*exp(-pi^2/ll^2*alpha*t); % Analytical function of temperature

dx = 0.001; % m
xx = (0:dx:ll)';
Nx = length(xx);

dtvett = [0.1 0.2 0.5 1 2 5 10 20];
errore = zeros(length(dtvett),1);

for zz = 1:length(dtvett)

    dt = dtvett(zz); % s
    time = (0:dt:tend);
    Nt = length(time);
    
    aa = alpha*dt/dx^2;
    
    Tm = T(xx);

    for ii = 2:Nt
    
        % Diagonals of the matrix at time m
        sub_diag_1 = -1/2*aa*ones(Nx,1);
        main_diag_1 = (1+aa)*ones(Nx,1);
        sup_diag_1 = sub_diag_1;
        
        AA_1 = spdiags([sub_diag_1 main_diag_1 sup_diag_1],-1:1,Nx,Nx);
        
        % Diagonals of the matrix at time m+1
        sub_diag_2 = 1/2*aa*ones(Nx,1);
        main_diag_2 = (1-aa)*ones(Nx,1);
        sup_diag_2 = sub_diag_2;
        
        AA_2 = spdiags([sub_diag_2 main_diag_2 sup_diag_2],-1:1,Nx,Nx);
        
        AA_1(1,1) = 1;
        AA_1(1,2) = 0;
        
        AA_1(end,end-1) = 0;
        AA_1(end,end) = 1;
        
        
        AA_2(1,1) = 1;
        AA_2(1,2) = 0;
        
        AA_2(end,end-1) = 0;
        AA_2(end,end) = 1;
    
        
        bb(1) = T0;
        bb(end) = T0;
        
        
        TT = AA_1\(AA_2*Tm);
        
        Tm = TT;
    
    end

Tanal = T_an(xx,tend); % analytical temperature
err = norm(TT-Tanal)/norm(Tanal-T0);
errore(zz) = err;

end

% figure(1)
% plot(xx,TT)
% uncomment to see the plot at t = 100 s

figure(2)
loglog(dtvett,errore,'-*k')
title('Relative error at different dt')
xlabel('dt size [s]')
ylabel('Relative error [-]')
legend('dx = 0.001 m','Location','northwest')
grid on

