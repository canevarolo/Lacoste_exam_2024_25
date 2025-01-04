%% VEDERE FUNZIONI SE SONO GIUSTE, ERRORE NEL GRAFICO RADIALE



% To study an innovative industrial process, a solid object made of a homogeneous material and shaped like a cylinder is constructed and tested in an
% induction furnace.
% The cylinder is, to a first approximation, cooled on all surfaces by an inert gas at a constant temperature of 20 °C with a heat transfer coefficient
% of 10 W/(m² K) (the anchoring system is neglected).  
% During nominal operation, the inductive heating deposits a total power of ~10 kW, distributed according to the function:  

% **qfunz = @(r,z) aa*(1-(r./RR).^2).*cos(pi/HH.*z)**  

% where *r* is the radial coordinate, *z* is the axial coordinate originating at the midpoint of the cylinder’s axis, the constant *aa* is 23 kW/m³,
% *R* is the radius (60 cm), and *H* is the cylinder’s height (120 cm).  

% The goal is to study the steady-state temperature distribution in the cylinder.
% To achieve this, a detailed and comprehensive report should be produced that includes:  

% a. The temperature distribution calculated using two 1D models (axial and radial, respectively).  
% b. A convergence study, including verification of the order of convergence, for one of the 1D cases treated in point a, as chosen.  



% Exam 31 january 2020
% Simone Canevarolo
% S269893
% 03/01/2025

clear all
close all
clc

HH = 120e-2; % m
RR = 60e-2; % m
Tgas = 293; % K
hh = 10; % W/m^2/K
aa = 23e3; % W/m^3
kk = 7; % W/m/K
rovol = 1600; % kg/m^3
cp = 2400; % J/kg/K

qtot = 10e3; % W

qfunz = @(r,z) aa*(1-(r./RR).^2).*cos(pi/HH.*z);

dz = 1e-4; % m
zz = (0:dz:HH/2)';
Nz = length(zz);

As = 2*pi*RR*HH/2; % m^2
Vol = pi*RR^2*HH/2; % m^3

% Axial solution

qfunz_axial = @(z) aa.*cos(pi/HH.*z);

sub_diag = kk/dz^2*ones(Nz,1);
main_diag = (-2*kk/dz^2-hh*As/Vol)*ones(Nz,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nz,Nz);

ii = 0;
zspace = 0;
bb = zeros(Nz,1);

for ii = 1:Nz
    
    zspace = ii*dz;
    bb(ii) = -qfunz_axial(zspace)-Tgas*hh*As/Vol;

end

AA(1,1) = -1;
AA(1,2) = 1;
bb(1) = 0;

AA(end,end-1) = -1;
AA(end,end) = 1+hh*dz/kk;
bb(end) = Tgas*hh*dz/kk;

TT = AA\bb;

figure(1)
plot(zz,TT_ax,'b','LineWidth',2)
title('Temperature along z axis')
xlabel('height [m]')
ylabel('Temperature [K]')


% Radius solution

dr = 1e-3; % m
rr = (0:dr:RR)';
Nr = length(rr);

As_rad = 2*pi*RR^2;

sub_diag = kk.*(1-dr./rr/2);
main_diag = (-2*kk-hh*As_rad/Vol)*ones(Nr,1);
sup_diag = kk.*(1+dr./rr/2);

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nr,Nr);

qfunz_rad = @(r) aa.*(1-(r/RR)^2);

jj = 0;
radspace = 0;
bb = zeros(Nr,1);

for jj = 1:Nr

    radspace = dr*jj;
    bb(jj) = -qfunz_rad(radspace)-Tgas*hh*As_rad/Vol;

end

AA(1,1) = -1;
AA(1,2) = 1;
bb(1) = 0;

AA(end,end-1) = -1;
AA(end,end) = 1+hh*dr/kk;
bb(end) = Tgas*hh*dr/kk;

TT = AA\bb;

figure(2)
plot(rr,TT)
