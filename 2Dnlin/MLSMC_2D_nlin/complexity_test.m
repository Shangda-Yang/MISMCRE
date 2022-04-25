close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver2Dnlin'));

sigma  = 0.5;
 
Eps = [ 0.001 0.0025 0.005 0.01 0.025 ];
% observation points
x_data = [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75]; 

loadpath = fullfile(filepart, 'Results','observations.mat');
load(loadpath,'data')

k = 2;
tic;
mlsmc_complexity(@opre_MLSMC_l,Eps,sigma,data,x_data,k);
toc;

rmpath(fullfile(filepart,'solver2Dnlin'));