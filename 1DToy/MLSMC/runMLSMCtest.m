% Test routine of MLSMC with
% the self-normalised increments estimator and ratio estimator,
% also computes the results for SLSMC for comparsion.

close all; clear all;
tic

[filepart,~,~] = fileparts(pwd); 
addpath(fullfile(filepart, 'pde_solver_1D'))

randn('seed',20)
rand('seed',20)

fprintf(1,'\n ---- Toy Example Ratio Estimator MLSMC ---- \n');
N      = 10000;       % samples for convergence tests
L      = 6;           % levels for convergence tests
Eps    = [ 0.0025 0.005 0.01 0.02 0.04 ];

% generate data
x_data = linspace(0, 1, 10); % corresponding point of data
sigma  = 0.2;
data   = -50*rand(1)*(x_data.^2-x_data) + sigma*randn(length(x_data),1);

filename = 'Toy_Example_MLSMC';
fp = fopen([filename '.txt'],'w');
mlsmc_test(@opre_mlsmc_l, N,L,Eps, fp,sigma, data, x_data);
fclose(fp);

toc

% plot results
% mlmc_plot(filename);


