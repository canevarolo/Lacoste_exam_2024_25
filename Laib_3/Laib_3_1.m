% In an infinite wall of thickness δ = 50 mm (see Fig. 1) with thermal conductivity k = 5 W/(m·K),
% there is a volumetric heat generation of [q''' = 500 kW/m³. Surface A is adiabatic, while surface B
% is maintained at a constant temperature T₀ = 300 K.
% Assuming the problem is 1D along x, analytically calculate and plot the temperature distribution within the wall.
% Then, recalculate the temperature distribution using the finite difference method and compare the numerical solution
% with the analytical one. Finally, verify energy conservation.

% Laib 3, exercise 1
% Simone Canevarolo
% S269893
% 28/10/2024

clear all
close all
clc

ll = 50e-3; % m
kk = 5; % W/m/K
qv = 500e3; % W/m^3
T0 = 300; % K

dx = 1e-5;
xx = (0:dx:ll)';
Nx = length(xx);

sub_diag = ones(Nx,1);
main_diag = -2*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
AA = spdiags(Band,-1:1,Nx,Nx);

bb = -qv*dx^2/kk*ones(Nx,1);

AA(1,1) = -1;
AA(1,2) = 1;
bb(1) = 0;

AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = T0;

TT = AA\bb;

plot(xx,TT,"r",'Linewidth',2)
hold on

AA(1,1) = -1;
AA(1,2) = 1;
bb(1) = qv*dx^2/2/kk;

plot(xx,TT,"k",'LineWidth',2)
hold off

% Energy conservation

flow_in = qv*ll;
flow_out = abs(kk/dx*(TT(end-1)-TT(end)));
err1 = abs(flow_out - flow_in)/abs(flow_in);
fprintf('Relative mistake is %.5f with method 1\n',err1)
flow_out_2 = abs(kk/dx*(TT(end-1)-TT(end))+qv*dx/2);
err2 = abs(flow_out_2 - flow_in)/abs(flow_in);
fprintf('Relative mistake is %.5f with method 2\n',err2)
