% Test de ode_integration
clear;
clear all;
[V,C] = ode_integration(1700, [-pi/8,0.1,0.1,0.1], 1500e3)