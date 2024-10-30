% In the same case of the exercise 1 (an infinite wall with 50 mm of thinkness and heat conducibility
% of 500 W/m/K and volumetric heat generation of 500 kW/m^3), the insulation is removed from the adiabatic surface,
% so that surface A is in contact with a fluid at temperature 273 K with which it exchanges
% heat by convection (heat transfer coefficient h_f = 100 W/(m² K)), see Fig. 4.
% Calculate analytically and numerically the temperature distribution in the wall,
% verify the incoming and outgoing heat fluxes, and compare the two solutions on the same graph.

% Laib 3, exercise 4
% Simone Canevarolo
% S269893
% 30/10/2024

clear all
close all
clc

ll = 50e-3; % m
kk = 5; % W/m/K
qv = 500e3; % W/m^3
T0 = 300; % K
hh = 100; % W/m^2/K
Tf = 273; % K

dx = 1e-4;
xx = (0:dx:ll)';
Nx = length(xx);

sub_diag = ones(Nx,1);
main_diag = -2*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = -qv*dx^2/kk*ones(Nx,1);

AA(1,1) = hh+kk/dx;
AA(1,2) = -kk/dx;
bb(1) = Tf*hh;

AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = T0;

TT = AA\bb;

figure(1)
plot(xx,TT,'r','LineWidth',2)
title('Temperature distribution')
xlabel('Length [m]')
ylabel('Temperature [K]')

% Energy conservation

flow_in = qv*ll;
flow_out = abs(hh*(Tf-TT(1)))+abs((T0-TT(end-1))/dx*kk);

err1 =abs(flow_out - flow_in)/abs(flow_in);

fprintf('The error is %.5f \n',err1)

figure(2)
flowbar = bar([flow_in,flow_out])
colors = [0 1 0; 1 0 0];
flowbar.FaceColor = 'flat';
flowbar.CData = colors;
title('Flow comparison');
ylabel('Flow [W]');
set(gca,'xticklabels',{'Flow in' 'Flow out'});