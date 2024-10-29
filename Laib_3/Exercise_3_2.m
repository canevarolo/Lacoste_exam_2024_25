% In an industrial facility, an insulated pipe (insulation thermal conductivity k = 0.15 W/(m·K), with an
% inner diameter Dₙ = 20 mm and an outer diameter Dₒ = 50 mm (see Fig. 2)) carries water at a constant
% temperature Tw = 15 °C, with an internal heat transfer coefficient hₙ = 1000 W/m²K.
% In the environment where it is located, the outer surface of the pipe receives a uniform heat flux q'' = 2 kW/m².
% Calculate and display the temperature distribution (under steady-state conditions) through the insulation
% (along the "r" direction) and compare it with the analytical solution. Finally, verify energy conservation

% Laib 3, exercise 2
% Simone Canevarolo
% S269893
% 29/10/2024

clear all
close all
clc

kk = 0.15; % W/m/K
Din = 20e-3; % m
Dout = 50e-3; % m
Tw = 288; % K
hh = 1000; % W/m^2/K
qout = 2e3; %W/m^2

% radius = (Dout-Din)/2;

dr = 1e-4; % m
rr = (Din/2:dr:Dout/2)';
Nr = length(rr);

sup_diag = 1+dr/2./rr;
main_diag = -2*ones(Nr,1);
sub_diag = 1-dr/2./rr;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
AA = spdiags(Band,-1:1,Nr,Nr);

bb = zeros(Nr,1);

AA(1,1) = -1-hh*dr/kk;
AA(1,2) = 1;
bb(1) = -Tw*hh*dr/kk;

AA(end,end-1) = -1;
AA(end,end) = 1;
bb(end) = qout*dr/kk;

TT = AA\bb;

plot(rr,TT-273,'LineWidth',2)

% Energy conservation

flow_in = qout;
flow_out = abs(kk/dr*(TT(end-1)-TT(end)));

err1 = abs(flow_out - flow_in)/abs(flow_in);

fprintf('Errore relativo è %.5f \n', err1)
