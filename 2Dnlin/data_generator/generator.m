%% pre-settings
close all; clear all;

randn('seed',20)
rand('seed',20)

sigma  = 0.5;
eps = 1e-10;
% observation points
x_data = [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75];
% % generate data
u = 2*rand(2,1) - 1;

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver2Dnlin'));
%% generate data
tic;
data =  data_generator(x_data,u,sigma);
toc;
savepath = fullfile(filepart, 'Results','observations.mat');
save(savepath,'data','u')
%% approximate reference solution
addpath(fullfile(filepart,'MLSMC_2D_nlin'));
k = 2;
tic;
[solutions] = mlsmc_solution(@opre_MLSMC_l,eps,sigma,data,x_data,k);
toc;
savepath = fullfile(filepart, 'Results','solutions.mat');
save(savepath,'solution');
rmpath(fullfile(filepart,'MLSMC_2D_nlin'))
rmpath(fullfile(filepart,'solver2Dnlin'));
%% data generator
% ------------------------------------------------------ %
% function [data] = data_generator(x_data, u, sigma)
% generate data based on the model y = G(u) + v
% input:  x_data = observation points
%         u      = random inputs
%         sigma  = std deviation of noise
% output: data   = data from the model
% ------------------------------------------------------ %
function [data] = data_generator(x_data, u, sigma)
lx = 10;
ly = 10;

Gu = pde_solver_2D(lx,ly,x_data,u);

data = Gu + sigma*randn(length(x_data), 1);
end
%% solution generator
% --------------------------------------------------------------------- %
% function [solutions] = mlsmc_solution(mlsmc_l,eps,sigma,data,x_data,k)
% inputs: mlsmc_l      = function for level l estimator
%         eps          = required accuracy
%         sigma        = std deviation of the error of y - G
%         data         = observations (values of y) vector
%         x_data       = corresponding observation points
%         K            = starting refinement level
% ouputs: solutions(1) = reference by sn
%         solutions(2) = reference by re 
% --------------------------------------------------------------------- %
function [solutions] = mlsmc_solution(mlsmc_l,eps,sigma,data,x_data,K)
% loading emperical results first
loadpath = fullfile(filepart, 'Results','mlsmc','emperical','mlsmc_emperical.mat');
load(loadpath,'se1','cost')

solutions = zeros(1,2);

vc1 = 2*sum(sum(sqrt(se1.*cost))).*sqrt(se1./cost);

L = ceil(log2(sqrt(2)/eps/3)/2)-K;
M = ceil(vc1./eps/eps);

sumlsn  = 0;
sumlret = 0;
sumlreb = 0;

for l = 0:L
    N = M(l+1);
    [sums,~] = mlsmc_l(l,l,N,sigma,data,x_data,K);
    sumlsn  = sumlsn  + sums(1);
    sumlret = sumlret + sums(2);
    sumlreb = sumlreb + sums(3);
end
solutions(1) = sumlsn;
solutions(2) = sumlret/sumlreb;

end
