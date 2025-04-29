% In un impianto industriale è presente un connettore metallico in lega di rame (ro = 5000 kg/m^3, cp = 100 J/kg/K, k = 200 W/m/K) di lunghezza 0.2 m,
% con sezione esagonale di lato 1 cm.
% Il connettore attraversa una membrana, di spessore trascurabile, che separa due ambienti: il primo a temperatura 50 °C, mentre il secondo a 20 °C.
% A causa delle vibrazioni legate al funzionamento dell'impianto, il connettore va periodicamente in contatto termico perfetto con un serbatoio a 100 °C
% (sulla base di sinsitra).
% Il contatto termico dura 20 secondi e si ripete una volta al minuto.
% Calcolare l'andamento di temperatura esterna a contatto con la faccia di sx, la temperatura della faccia di sinistra, la temperatura al centro del connettore
% e la temperatura sulla faccia di dx. Da calcolore, inoltre, l'andamento delle medesime a stazionario.


% In an industrial plant, there is a metallic connector made of a copper alloy (ρ = 5000 kg/m³, cₚ = 100 J/kg·K, k = 200 W/m·K), with a length of 0.2 m and a
% hexagonal cross-section with a side of 1 cm.
% The connector passes through a membrane of negligible thickness that separates two environments: one at a temperature of 50 °C and the other at 20 °C.
% Due to vibrations related to the operation of the plant, the connector periodically comes into perfect thermal contact with a reservoir at 100 °C (on the
% left-hand side).
% The thermal contact lasts for 20 seconds and repeats once every minute.
% Calculate the temperature trend at the external contact on the left face, the temperature of the left face itself, the temperature at the center of the connector
% and the temperature on the right face. Also calculate the steady-state behavior of the same temperatures.




clear all
close all
clc

rovol = 5e3; % kg/m^3
cp = 100; % J/kg/K
kk = 200; % W/m/K

ll = 0.2; % m
lato = 1e-2; % m

T1 = 50+273; % K
T2 = 20+273; % K
Tcon = 100+273; % K

hh1 = 20; % W/m^2/K
hh2 = 40; % W/m^2/K

dx = 1e-3; % m
xx = (0:dx:ll)';
Nx = length(xx);

dt = 1; % s

aa = kk*dt/dx^2/rovol/cp;

% Tm = zeros(Nx,1); % K
% Tm(1:Nx/2) = T1; % K
% Tm(Nx/2+1:Nx) = T2; % K

Tm = Tcon*ones(Nx,1);

As = 6*lato*ll; % m^2, Thermal exchange area
Vol = 6*((sqrt(3)/2)*lato*lato/2)*ll; % m^3, Volume

hh = hh1*(xx<=ll/2)+hh2*(xx>ll/2); % W/m^2/K
Tair = T1*(xx<=ll/2)+T2*(xx>ll/2); % K

tper = 60; % s, time of the period
per = 1; % number of periods passed
ttot = 8*tper; % s
Tright = zeros(tper,1);
tempo = 0:dt:ttot;
time = 0;

sub_diag = -aa*ones(Nx,1);
main_diag = 1+2.*aa+hh.*As./Vol*dt/rovol/cp;
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = Tair.*hh*As/Vol*dt/rovol/cp;

AA(1,1) = -kk/dx;
AA(1,2) = kk/dx+hh1;
bb(1) = T1*hh1;

AA(end,end-1) = -kk/dx;
AA(end,end) = kk/dx+hh2;
bb(end) = T2*hh2;

TT = AA\bb;

Trightmax = T1;

for ii = 1:length(tempo)
   
    time = time+dt;

    sub_diag = -aa*ones(Nx,1);
    main_diag = 1+2.*aa+hh.*As./Vol*dt/rovol/cp;
    sup_diag = sub_diag;
    
    Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
    
    AA = spdiags(Band,-1:1,Nx,Nx);
    
    bb = Tm+Tair.*hh*As/Vol*dt/rovol/cp;

    if time >= 20 && time < 40

        % Dirichlet 

        AA(1,1) = 1;
        AA(1,2) = 0;
        bb(1) = Tcon;
    
        AA(end,end-1) = -kk/dx;
        AA(end,end) = kk/dx+hh2;
        bb(end) = T2*hh2;

    else

        % Robin

        AA(1,1) = -kk/dx;
        AA(1,2) = kk/dx+hh1;
        bb(1) = T1*hh1;
    
        AA(end,end-1) = -kk/dx;
        AA(end,end) = kk/dx+hh2;
        bb(end) = T2*hh2;

    end

    TT = AA\bb;

    Tleft(ii) = TT(1);
    Tmiddle(ii) = TT(round(Nx/2));
    Tright(ii) = TT(end);

    if time == tper

        time = 0;
        per = per+1;

    end

    if per > 5
        if Tright(ii) > Trightmax

            Trightmax = Tright(ii);
            Tsave = TT;

        end
    end

    Tm = TT;

end

figure(1)
plot(tempo,Tright-273,'b')
hold on
plot(tempo,Tleft-273,'r')
hold on
plot(tempo,Tmiddle-273,'g')

figure(2)
plot(xx,Tsave-273)





dxvett = [5e-2 1e-2 5e-3 1e-3 5e-4 1e-4];

errdx = zeros(length(zz),1);


for zz = 1:length(dxvett)





    errdx(zz) = 

end

figure(3)
loglog()



dtvett = [10 5 2 1 0.5 0.2 0.1 0.05];

errdt = 


for jj = 1:length(dtvett)

    dt = dtvett(jj);

    tper = 60; % s, time of the period
    per = 1; % number of periods passed

    ttot = 8*tper; % s
    Tright = zeros(tper,1);
    tempo = 0:dt:ttot;
    time = 0;

    for ii = 1:length(tempo)

        time = time+dt;

        sub_diag = -aa*ones(Nx,1);
        main_diag = 1+2.*aa+hh.*As./Vol*dt/rovol/cp;
        sup_diag = sub_diag;

        Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

        AA = spdiags(Band,-1:1,Nx,Nx);

        bb = Tm+Tair.*hh*As/Vol*dt/rovol/cp;

        if time >= 20 && time < 40

            % Dirichlet 

            AA(1,1) = 1;
            AA(1,2) = 0;
            bb(1) = Tcon;

            AA(end,end-1) = -kk/dx;
            AA(end,end) = kk/dx+hh2;
            bb(end) = T2*hh2;

        else

            % Robin

            AA(1,1) = -kk/dx;
            AA(1,2) = kk/dx+hh1;
            bb(1) = T1*hh1;

            AA(end,end-1) = -kk/dx;
            AA(end,end) = kk/dx+hh2;
            bb(end) = T2*hh2;

        end

        TT = AA\bb;

        Tright(ii) = TT(end);

        Tm = TT;

    end


    errdt(jj) = norm(Trightmax-Tright)/norm(Tright);

end

figure(4)
loglog(dtvett,errdt)
