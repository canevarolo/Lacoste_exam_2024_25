% Thermal energy storage systems often involve packed beds of solid spheres through which hot
% gas flows during the 'charging' phase and cold gas flows during the 'discharging' phase.
% In the charging process (shown in Fig. 2), heat transfer from the hot gas increases the
% thermal energy stored in the spheres, which are initially at a lower temperature; during the
% discharging phase, the stored thermal energy decreases as it is transferred from the
% hot spheres to the cooler gas.  
% Consider a bed composed of aluminum spheres with a diameter of 75 mm (density 2700 kg/m³,
% specific heat capacity 950 J/kg/K, thermal conductivity 240 W/m/K) and a charging process
% where the gas enters the storage unit at a temperature of 300 degrees Celsius.
% If the initial temperature of the spheres is 25 degrees Celsius and the heat transfer
% coefficient is 75 W/m²K, numerically determine, after verifying that the problem can be
% addressed with a lumped-parameter model, how long it takes for the system to store 90 percent
% of the maximum storable thermal energy.
% Finally, evaluate whether and how much time would be saved in storing the
% same amount of energy if the spheres were made of copper (density 8900 kg/m³,
% specific heat capacity 400 J/kg/K) instead of aluminum, while maintaining the same size.

% Laib 8, exercise 4
% Simone Canevarolo
% S269893
% 27/11/2024

clear all
close all
clc

dd = 75e-3; % m
rovol_all = 2700; % kg/m^3
cp_all = 950; % J/kg/K
kk_all = 240; % W/m/K
Tgas = 300+273; % K
Tin = 25+273; % k
hh = 75; % W/m^2/K

As = 4*pi*dd^2/4; % area for heat exchange
Vol = 4/3*pi*(dd/2)^3; % volume of the sphere

% Biot number (should be Bi<0.1 in order to use the zero-dimensional
% problem)

Bi = hh*(dd/2)/kk_all;

% Now, we can start to code since the condition is satisfied.
% We consider the max of energy that can be stored as if the aluminium
% balls could have the same temperature as the air flow.
% The code will be solved in Forward Euler (FE), explicit method

Emax_theory_all = @(T) rovol_all*Vol*cp_all*(T-Tin);
Emax = Emax_theory_all(Tgas);

TT = Tin;
Tm = Tin;
Energy = rovol_all*Vol*cp_all*(TT-Tin);

err = 1; % error
toll = 1e-5; % tollerance
dt = 1; % s
time = 0; % s
ii = 1; % iterations

TTvett_all(1) = TT;

aa = dt*As*hh/rovol_all/cp_all/Vol;

while Energy < 0.9*Emax && ii<1000

    ii = ii+1;
    time = time+dt;

    TT = Tm+aa*(-Tm+Tgas);

    TTvett_all(ii) = TT;

    Energy = rovol_all*Vol*cp_all*(TT-Tin);

    Tm = TT;

end

tt = (0:dt:time)';

figure(1)
plot(tt,TTvett_all-273,'k','LineWidth',2)
title('Temperature on alluminium vs time')
xlabel('time [s]')
ylabel('Temperature [°C]')
% legend(Alluminium temperature)

Time_all = time;


%% COPPER

% This time, the problem will be solved using Backward Euler (BE), implicit
% method

rovol_cu = 8900; % Kg/m^3
cp_cu = 400; % W/m/K
% same dimentions

Emax_theory_cu = @(T) rovol_all*Vol*cp_all*(T-Tin);
Emax = Emax_theory_cu(Tgas);

TT = Tin;
Tm = Tin;
Energy = rovol_cu*Vol*cp_cu*(TT-Tin);

time = 0; % s
ii = 1;

TTvett_cu(1) = TT;

aa = dt*As*hh/rovol_cu/cp_cu/Vol;

while Energy < 0.9*Emax && ii<1000

    ii = ii+1;
    time = time+dt;

    TT = (Tm+aa*Tgas)/(1+aa);

    TTvett_cu(ii) = TT;

    Energy = rovol_cu*Vol*cp_cu*(TT-Tin);

    Tm = TT;

end

tt = (0:dt:time)';

figure(2)
plot(tt,TTvett_cu-273,'r','LineWidth',2)
title('Temperature on copper vs time')
xlabel('time [s]')
ylabel('Temperature [°C]')
% legend(Copper temperature)


timesave = abs(time-Time_all);
fprintf('\nThe time saved is %.5f seconds\n', timesave)
