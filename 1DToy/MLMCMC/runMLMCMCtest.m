% Test routine for MLMCMC

close all; clear all;

[filepart,~,~] = fileparts(pwd); 
addpath(fullfile(filepart, 'pde_solver_1D'))

randn('seed',20)
rand('seed',20)

fprintf(1,'\n ---- Toy Example MLMCMC ---- \n');
N      = 10000;       % samples for convergence tests
L      = 6;           % levels for convergence tests
Eps    = 2.^(-11:1:-4);

% generate data
x_data = linspace(0, 1, 10); % corresponding point of data
sigma  = 0.2;
data   = -0.5*rand(1)*(x_data.^2-x_data) + sigma*randn(length(x_data),1);

filename = 'Toy_Example_MLMCMC';
fp = fopen([filename '.txt'],'w');
mlmcmc_test(@opre_mlmcmc_l, N,L,Eps, fp,sigma, data, x_data);
fclose(fp);


% plot results
% mlmc_plot(filename);
