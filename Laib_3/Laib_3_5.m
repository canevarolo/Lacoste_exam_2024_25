% Pressurized water at a constant temperature of Tw = 400 K flows through 
% a cylindrical pipe (inner diameter Din = 20 mm, outer diameter Dout = 24 mm,
% conductivity k = 0.035 W/(m K), heat transfer coefficient hin = 100 W/(m² K)).
% On the outer surface, the pipe is cooled by air at Ta = 300 K (heat transfer coefficient hout = 10 W/(m² K)).
% Calculate analytically and numerically the radial temperature profile across the cross-section
% of the pipe and compare the two results in a graph.
% Verify the conservation of energy.

% Laib 3, exercise 5
% Simone Canevarolo
% S269893
% 01/11/2024

clear all
close all
clc

Tw = 400; % W
Din = 20e-3; % m
Dout = 24e-3; % m
kk = 0.035; % W/m/K
hin = 100; % W/m^2/K
Ta = 300; % K
hout = 10; % W/m^2/K

dr = 1e-5;
rr = (Din/2:dr:Dout/2)';
Nr = length(rr);

sup_diag = 1+dr./2./rr;
main_diag = -2*ones(Nr,1);
sub_diag = 1-dr./2./rr;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nr,Nr);

bb = zeros(Nr,1);

AA(1,1) = hin+kk/dr;
AA(1,2) = -kk/dr;
bb(1) = hin*Tw;

AA(end,end-1) = -1;
AA(end,end) = 1+hout*dr/kk;
bb(end) = Ta*hout*dr/kk;

TT = AA\bb;

figure(1)
plot(rr*1e3,TT,'r','LineWidth',2)
title('Temperature along the radius')
xlabel('Thickness [mm]')
ylabel('Temperature [K]')

% Energy conservation

flow_in =  abs(hout*(TT(end)-Ta))*pi*Dout;
flow_out = abs(hin*(TT(1)-Tw))*pi*Din;

err1 = abs(flow_out- flow_in)/abs(flow_in);

fprintf('The error is %.5f', err1)