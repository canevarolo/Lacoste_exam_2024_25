% Datas from Laib 9_1_a
% 
% Recalculate the temperature distribution in the rod using finite differences
% with an explicit time-stepping scheme. Compare the solutions obtained at
% t = 100 s for Δx = 0.5 cm and Δt = 1.0 s, 1.25 s, and 2 s in Figure #4.
% What happens to the solution as Δt increases? Why?

% Laib 9, exercise 1 D
% Simone Canevarolo
% 11/12/2024
% S269893

clear all
close all
clc

ll = 10e-2; % length, m
alpha = 1e-5; % thermal diffusivity, m^2/s
tzero = 0; % time to begin, s

T0 = 300; % temperature of both borders, K
tend = 100; % time to end, s

T = @(x) T0+50*sin(pi*x/ll); % temperature function

T_an = @(x,t) T0+50*sin(pi*x/ll)*exp(-pi^2/ll^2*alpha*t); % Analytical function of temperature

dx = 0.005; % m
xx = (0:dx:ll)';
Nx = length(xx);

figure(1)
plot(xx,T_an(xx,tend));
hold on

Tm = T(xx);

dtvett = [1, 1.25, 2];

for zz = 1:length(dtvett)

dt = dtvett(zz);
time = (0:dt:tend);
Nt = length(time);

    for ii = 1:Nt

        aa = alpha*dt/dx^2;

        sub_diag = aa*ones(Nx,1);
        main_diag = (1-2*aa)*ones(Nx,1);
        sup_diag = sub_diag;
        
        Band = [sub_diag main_diag sup_diag];
        
        AA = spdiags(Band,-1:1,Nx,Nx);
                
        bb = Tm;
    
        AA(1,1) = 1;
        AA(1,2) = 0;
        bb(1) = T0;
    
        AA(end,end-1) = 0;
        AA(end,end) = 1;
        bb(end) = T0;
    
        TT = AA*bb;
    
        Tm = TT;
    
    end

plot(xx,TT)
hold on

end

title('Temperature with explicit method')
xlabel('Lenght [m]')
ylabel('Temperature [K]')
legend('Analitic','dt = 1 s', 'dt = 1.25 s','dt = 2 s','Location','best')

% As can be observed, the system is not well-conditioned, and for Δt = 2 s,
% we fail to satisfy the absolute stability condition (α · Δt / Δx² ≤ 0.5).
% A graph like the one visible at that time is clearly not representative of reality.  
% The Forward Euler method is, in fact, conditionally absolutely stable.