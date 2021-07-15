%{ 
    CU Sounding Rocket Lab - Reaper Engine Development
    Thrust Chamber Thermostructural Analysis
    
    Samuel Ciesielski
%}

clear; close all; clc;

import @engine

%% SETUP

%%% NOZZLE PARAMETERS (GRCop-42)
n = 200;                                                                    % number of 

t_w         = ;                                                              % [m] nozzle wall thickness
lambda_w    = 344;                                                           % [W/(m*K)] thermal conductivity 
E           = ;                                                              % [] modulus of elasticity
a           = ;                                                              % thermal expansion coefficient
v           = ;                                                              % Poisson's ratio

nozParams = [t_w lambda_w E a v];



%%% ENGINE PERFORMANCE PARAMETERS
p_c         = ;                                                             % [N/m^2] chamber pressure
T_0         = ;                                                             % [K] stagnation temp
cp_g        = ;                                                             % [J/(K*kg)] specific heat of core flow at constant pressure
gamma       = ;                                                             % specific heat ratio
wDot        = ;                                                             % [N/s] total prop weigth flowrate

perfParams = [p_c, T_0, cp_g, gamma, wDot];



%%% COOLING SYSTEMS PARAMETERS
n_cc        = ;                                                             % number of cooling channels
T_c_i       = ;                                                             % [K] regen coolant input temp
dm_c        = ;                                                             % [kg/s] fuel/regen coolant massflow rate
h_cc        = [];                                                           % [m] channel height at [inlet throat injector]
w_cc        = [];                                                           % [m] channel inlet  at [inlet throat injector]
t_j         = ;                                                             % [m] jacket thickness
t_f         = ;                                                             % [m]fin thickness

dm_f        = ;                                                             % [kg/s] film massflow
lInj        = ;                                                             % [m] axial position of film injection

lambda_c    = ;
Q_cVap      = ;
cp_c        = ;
mu_c        = ;
k_c         = ;

cProperties = [lambda_c, Q_cVap, cp_c, mu_c, k_c];

%%% CONSTRUCT MODEL

Reaper      = engine(perfParams, nozParams, "rpa", "Engine2_Thermal.txt", n);
Reaper      = applyRegen(Reaper, T_c_i, dm_c, h, w, [t_f t_j], cProperties, n_cc);
Reaper      = applyFilm(Reaper, dm_f, lInj);
%Reaper     = applyTBC(Reaper, t_tbc, lambda_tbc);


%% STEADY STATE NOZZLE TEMPS & STRESSES

Reaper      = Bartz(Reaper);
Reaper      = nozStresses(Reaper);

plotThermal(Reaper);
%plotStresses(Reaper);

