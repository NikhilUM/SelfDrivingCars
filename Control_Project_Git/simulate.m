clc
clear
close all
% forward Integrates Simulation
[Y_sim,U_sim,t_total,t_update,Xobs] = forwardIntegrate();
save('refs/sim_result2.mat');

