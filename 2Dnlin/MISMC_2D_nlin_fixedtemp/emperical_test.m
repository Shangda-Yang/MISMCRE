close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver2Dnlin'));

N      = 1000;        % samples for convergence tests
Lx     = 4;           % levels for convergence tests
Ly     = 4;
sigma  = 0.5;

% observation points
x_data = [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75]; 

loadpath = fullfile(filepart, 'Results','observations.mat');
load(loadpath,'data')

k = 2;
tic;
mismc_emperical(@opre_MISMC_l,N, Lx, Ly,sigma,data,x_data,k);
toc;