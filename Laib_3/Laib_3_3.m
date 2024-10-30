% A cylindrical copper rod (diameter D = 0.01 m, length L = 4 m) carries a current I = 1 kA.
% The rod (see Fig. 3) is immersed in (and cooled by) a liquid nitrogen bath at Tb = 77 K
% (heat transfer coefficient h = 500 W/(m² K)).
% Both ends are maintained at the constant temperature T = Tb. Numerically calculate the
% steady-state temperature distribution along the rod in a script, verify energy conservation,
% and plot the temperature along the rod in a figure. Perform the exercise considering the minimum
% computational domain.
% (Copper electrical resistivity ρel = 1.75e-8 Ω m, Copper thermal conductivity k = 350 W/(m K)

% Laib 3, ecercise 3
% Simone Canevarolo
% S269893
% 30/10/2024


% FULL DOMAIN CASE

clear all
close all
clc

ll = 4; % m
DD = 1e-2; % m
II = 1e3; % A
Tb = 77; % K
hh = 500; % W/m^2/K
roel = 1.75e-8; % Ohm*m
kk = 350; % W/m/K

Abase = pi*DD^2/4; % area of the base
As = pi*DD*ll; % area of the lateral surface of exchange
Res = roel*ll/Abase; % electrical resistance
Vol = Abase*ll; % volume of the cable
qv = II^2*Res/Vol; % heat generation out of volume

dx = 5e-3; % m
xx = (0:dx:ll)';
Nx = length(xx);

cc = hh*As/Vol*dx^2/kk; % term of moltiplication

sub_diag = ones(Nx,1);
main_diag = (-2-cc)*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = (-qv*dx^2/kk-Tb*cc)*ones(Nx,1);

AA(1,1) = 1;
AA(1,2) = 0;
bb(1) = Tb;

AA(end,end-1) = -1;
AA(end,end) = 1;
bb(end) = 0;

TT = AA\bb;

figure(1)
plot(xx,TT,'r','LineWidth',2)
title('Temperature distribution along the full domain')
xlabel('Length [m]')
ylabel('Temperature [K]')

% Energy conservation

flow_in = qv*Abase*ll;
flow_out = kk/dx*(TT(1)-TT(2))*Abase+hh*pi*DD*trapz(xx,TT-Tb);
err1 = abs(flow_out - flow_in)/abs(flow_in);

fprintf('The relative error is %.5f', err1)

figure(2)
flowbar = bar([flow_in,flow_out])
colors = [0 1 0; 1 0 0];
flowbar.FaceColor = 'flat';
flowbar.CData = colors;
title('Flow comparison');
ylabel('Flow [W]');
set(gca,'xticklabels',{'Flow in' 'Flow out'});


%% HALF DOMAIN CASE


clear all
close all
clc

ll = 4; % m
DD = 1e-2; % m
II = 1e3; % A
Tb = 77; % K
hh = 500; % W/m^2/K
roel = 1.75e-8; % Ohm*m
kk = 350; % W/m/K
llhalf = ll/2;

Abase = pi*DD^2/4; % area of the base
As = pi*DD*llhalf; % area of the lateral surface of exchange
Res = roel*llhalf/Abase; % electrical resistance
Vol = Abase*llhalf; % volume of the cable
qv = II^2*Res/Vol; % heat generation out of volume

dx = 5e-3; % m
xx = (0:dx:llhalf)';
Nx = length(xx);

cc = hh*As/Vol*dx^2/kk; % term of moltiplication

sub_diag = ones(Nx,1);
main_diag = (-2-cc)*ones(Nx,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = (-qv*dx^2/kk-Tb*cc)*ones(Nx,1);

AA(1,1) = 1;
AA(1,2) = 0;
bb(1) = Tb;

AA(end,end-1) = -1;
AA(end,end) = 1;
bb(end) = 0;

TT = AA\bb;

figure(3)
plot(xx,TT,'k','LineWidth',2)
title('Temperature distribution along half of domain')
xlabel('Length [m]')
ylabel('Temperature [K]')

% Energy conservation

flow_in = qv*Abase*llhalf;
flow_out = kk/dx*(TT(1)-TT(2))*Abase+hh*pi*DD*trapz(xx,TT-Tb)
err1 = abs(flow_out - flow_in)/abs(flow_in);

fprintf('The relative error is %.5f \n', err1)

figure(4)
flowbar = bar([flow_in,flow_out])
colors = [0 1 0; 1 0 0];
flowbar.FaceColor = 'flat';
flowbar.CData = colors;
title('Flow comparison');
ylabel('Flow [W]');
set(gca,'xticklabels',{'Flow in' 'Flow out'});