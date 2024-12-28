% Steel spheres with a diameter of 12 mm (conductivity = 40 W/m/K, density = 7800 kg/m³, and specific heat = 600 J/kg/K)
% are tempered by rapidly heating them to 1150 K and then cooling them slowly to a final temperature of 400 K
% in an environment with air, whose temperature T_air increases over time: T_air = 325 K + 0.0375 K/s × t,
% where t is the time elapsed since the beginning of the cooling process.
% Calculate numerically, using the implicit Euler scheme, the time evolution of the temperature of the spheres,
% assuming a heat transfer coefficient of 20 W/m²K and that the entire external surface of the spheres is in
% contact with the air.
% Perform a one-dimensional analysis to show the time evolution of the average temperature of the spheres and
% determine after how long it becomes equal to that of the air. Verify the obtained result using the forward Euler scheme.
% Finally, compare the obtained results with the zero-dimensional analysis conducted in Exercise 4 of Lab #7.

% Laib 9, exercise 5
% Simone Canevarolo
% S269893
% 28/12/2024

clear all
close all
clc

% TERMINARE
