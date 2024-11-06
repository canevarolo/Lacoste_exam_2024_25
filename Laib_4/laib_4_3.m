% The cross-section of a copper cable (thermal conductivity kCu = 150 W/(m K),
% specific heat cp,Cu = 350 J/(kg K), density ρCu = 8900 kg/m^3, electrical resistivity
% ρel,Cu = 2 × 10^-8 Ω m) used for the combined transport of liquid nitrogen and electric
% power is shown in the figure. The conductor carries a current I = 20 kA and is internally
% cooled by the flow of nitrogen, which is at a temperature TN = 77 K (heat transfer
% coefficient hint = 1200 W/m^2K), while the outer side is thermally insulated by a layer
% of thermal insulation (thermal conductivity kins = 0.5 W/mK, specific heat cp,ins = 1800 J/kgK,
% density ρins = 2500 kg/m^3, with assumed infinite electrical resistivity), which, in turn,
% exchanges heat with stagnant air at Ta = 25 °C (heat transfer coefficient hext = 15 W/m^2K).
% The contact resistance between the internal copper cylinder and the insulation can be considered negligible.
% Calculate numerically and graphically represent the radial temperature distribution across the
% cross-section of the cable (both in the copper and in the insulation).
% By parametrically varying the thickness of the insulation between 0.1 mm and 15 mm, numerically calculate
% and graphically represent the value of the heat transferred per unit length through the insulation,
% with the aim of evaluating the critical insulation thickness t_c that minimizes the heat losses of the cable.
% Remember that the analytical expression for the critical radius is t_c = k_ins / h_ext - r_in, where r_in is
% the radius of the thermal insulation.


% Laib 4, exercise 3
% Simone Canevarolo
% S269893
% 5/11/2024

clear all
close all
clc

rint = 15e-3; % m, internal radius with N2
sscu = 10e-3; % m, thickness of copper
ssins = 1e-3; % m, thickness of insulant layer
rtot = rint+sscu+ssins; % m
sstot = sscu+ssins; % m

% Copper
kcu = 150; % W/m/L, heat conduction
cpcu = 350; % J/Kg/K, specific heat
rocu = 8900; % Kg/m^3, density
roelcu = 2e-8; % Ohm*m, electrical resistance

II = 20e3; % A

Tn = 77; % K, temperature of nitrogen
hint = 1200; % W/m^2/K

% insulat layer
kins = 0.5; % W/m/K
cpins = 1800; % J/Kg/K
roins = 2500; % Kg/m^3
% roelins = inf. for Hp.

Ta = 298; % K, air temperature externally
ha = 15; % W/m^2/K

dr = 1e-5;
rr = (rint:dr:rtot)';
Nr = length(rr);

rrcu = (rint:dr:(rint+sscu))';
Nrcu = length(rrcu);

Abasecu = pi*((rint+sscu)^2-rint^2);
qv = II^2*roelcu/Abasecu^2;

sub_diag = 1-dr./(2*rr);
main_diag = -2*ones(Nr,1);
sup_diag = 1+dr./(2*rr);

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nr,Nr);

bbcu = -qv.*dr^2./kcu*ones(Nrcu,1);
bbis = zeros(Nr-Nrcu,1);
bb = [bbcu; bbis];

AA(1,1) = 1+hint*dr/kcu;
AA(1,2) = -1;
bb(1) = Tn*hint*dr/kcu+qv*dr^2/2;

AA(Nrcu,Nrcu-1) = kcu;
AA(Nrcu,Nrcu) = -kcu-kins;
AA(Nrcu,Nrcu+1) = kins;
bb(Nrcu) = 0;

AA(end,end-1) = -1 ;
AA(end,end) = 1+ha*dr/kins;
bb(end) = Ta*ha*dr/kins;

TT = AA\bb;

figure(1)
plot(rr*1e3,TT,'r','LineWidth',2)
hold on
xline(25,'k','LineWidth',2)
hold off
xlim([15 26])
title('Temperature along the layers')
xlabel('Thickness [mm]')
ylabel('Temperature [K]')
legend('Temperature','Interface','Location','best')

%%

tc_1 = kins/ha-(rint+sscu);
qloss_1 = 2*pi*(rint+sscu+ssins)*ha*(TT(end)-Ta);

ssis = linspace(1e-4,15e-3,30);
qlosses = zeros(length(ssis),1);
tc = zeros(length(ssis),1);

for jj = 1:length(ssis)
    
    tc(jj) = kins/ha-ssis(jj); % critical thickness

    rtot = rint+sscu+ssis(jj); % m
    sstot = sscu+ssins; % m

    dr = 1e-5;
    rr = (rint:dr:rtot)';
    Nr = length(rr);
    
    rrcu = (rint:dr:(rint+sscu))';
    Nrcu = length(rrcu);
    
    Abasecu = pi*((rint+sscu)^2-rint^2);
    qv = II^2*roelcu/Abasecu^2;
    
    sub_diag = 1-dr./(2*rr);
    main_diag = -2*ones(Nr,1);
    sup_diag = 1+dr./(2*rr);
    
    Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
    
    AA = spdiags(Band,-1:1,Nr,Nr);
    
    bbcu = -qv.*dr^2./kcu*ones(Nrcu,1);
    bbis = zeros(Nr-Nrcu,1);
    bb = [bbcu; bbis];
    
    AA(1,1) = 1+hint*dr/kcu;
    AA(1,2) = -1;
    bb(1) = Tn*hint*dr/kcu+qv*dr^2/2;
    
    AA(Nrcu,Nrcu-1) = kcu;
    AA(Nrcu,Nrcu) = -kcu-kins;
    AA(Nrcu,Nrcu+1) = kins;
    bb(Nrcu) = 0;
    
    AA(end,end-1) = -1 ;
    AA(end,end) = 1+ha*dr/kins;
    bb(end) = Ta*ha*dr/kins;
    
    TT = AA\bb;

    qlosses(jj)=2*pi*(rint+sscu+ssis(jj))*ha*(TT(end)-Ta);

end

figure(2)
plot(ssis*1e3,qlosses,'-r*','LineWidth',2)
title('Heat loss along the diameter vs increase thickness insulant')
xlabel('Thickness [mm]')
ylabel('Heat loss [W/m]')
hold on
xline(tc_1*1e3,'LineWidth',2)
hold off
legend('Heat loss','Location','best')
