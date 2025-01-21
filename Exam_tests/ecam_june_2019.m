% A nuclear reactor fuel rod (L ~ 1 m) is composed of two layers, the fuel (MOX) and the cladding,
% in perfect thermal contact, as shown in Fig. 1.
% The rod, under nominal conditions, generates a uniform power of 3.5 MW/m³ (only in the fuel region)
% and is cooled by liquid sodium (Na) at 527.6 °C, with a heat transfer coefficient of 2000 W/(m²·K).
% Starting from the nominal steady-state condition, a transient begins in the reactor, inducing a
% power ramp according to the law:  
%
% qvol(t) = qvol,0 + b·t, where b = 310 kW/(m³·s).  
% 
% The power ramp ends at time ts, when the temperature at the center of the rod reaches 700 °C.
% After ts, the reactor is shut down, and power is generated by decay heat in both the fuel and the cladding,
% according to a given decay heat law:
% 
% qvolspento = @(t) qvolmax*((t-tend).^(-0.2)-t.^(-0.2))
% 
% Prepare a detailed report on the thermal transient analysis of the fuel rod, including:
% 
% a. The temperature distribution at the start of the transient and at time ts.  
% b. A grid independence study of the results for the nominal steady-state condition.  
% c. The evolution of the peak temperature in the fuel and cladding (using a logarithmic scale for the time axis)
% until reaching a new steady state.  
% d. A time discretization independence study of the results.  

% Simone Canevarolo
% S269893
% Exam test - June 2019

clear all
close all
clc


din = 8e-3; % m
dout = 10e-3; % m

raggio = dout/2; % m

rcomb = din/2; % m
rguaina = (dout-din)/2; % m

Tna = 527.6+273; % K
Tend = 700+273; % K
hh = 2000; % W/m^2/K
qvolzero = 3.5e6; % W/m^3

b = 310e3; % W/m^3/s

qvol_in = @(t) qvolzero+b.*t;

qvolspento = @(t) qvolmax*((t-tend).^(-0.2)-t.^(-0.2));

kkmox = 4.03; % W/m/K
romox = 11027; % kg/m^3
cpmox = 237; % J/kg/K

kkguaina = 16.3; % W/m/K
roguaina = 6753; % kg/m^3
cpguaina = 327; % J/kg/K

dr = 1e-5; % m
rr = (0:dr:raggio)';
rr1 = (0:dr:rcomb)';
rr2 = (rcomb:dr:raggio)';

Nr = length(rr);
Nr1 = length(rr1);
Nr2 = length(rr2);

dt = 0.1; % s
time = 0; % s

Tm = Tna*ones(Nr,1);
Tm1 = Tna*ones(Nr1,1);
Tm2 = Tna*ones(Nr2,1);

Tmax = Tna;

while Tmax < Tend

    aa1 = kkmox*dt/romox/cpmox/dr^2;
    aa2 = kkguaina*dt/roguaina/cpguaina/dr^2;
    
    sub_diag1 = -aa1.*(1-dr/2./rr1);
    main_diag1 = (1+2*aa1)*ones(Nr1,1);
    sup_diag1 = -aa1.*(1+dr/2./rr1);
    
    sub_diag2 = -aa2.*(1-dr/2./rr2);
    main_diag2 = (1+2*aa2)*ones(Nr2,1);
    sup_diag2 = -aa2.*(1+dr/2./rr2);
    
    sub_diag = [sub_diag1; sub_diag2];
    main_diag = [main_diag1; main_diag2];
    sup_diag = [sup_diag1; sup_diag2];
    
    
    Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
    
    AA = spdiags(Band,-1:1,Nr,Nr);
    
    bb = zeros(Nr,1);
    bb(1:Nr1) = Tm(1:Nr1)+qvol_in(time)*dt/romox/cpmox;
    bb(Nr1+1:Nr) = Tm(Nr1+1:Nr);

    AA(1,1) = -1;
    AA(1,2) = 1;
    bb(1) = 0;

    AA(Nr1,Nr1-1) = kkguaina;
    AA(Nr1,Nr1) = -kkguaina-kkmox;
    AA(Nr1,Nr1+1) = kkmox;
    bb(Nr1) = 0;

    AA(end,end-1) = -kkguaina/dr;
    AA(end,end) = kkguaina/dr+hh;
    bb(end) = Tna*hh;
    
    TT = AA\bb;
    
    Tmax = max(TT);
    
    time = time+dt;
    Tm = TT;

end

figure(1)
plot(rr,TT-273,'LineWidth',2)
