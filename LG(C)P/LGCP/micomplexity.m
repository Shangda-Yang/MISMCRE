close all; clear all;

randn('seed',20)
rand('seed',20)

% s2:   \sigma^2
% mu:   \theta_1
% a:    \theta_2
% b:    \theta_3
% rate: \beta^{\prime}

s2=1.91;
b=(33/pi)^2;
% a=s2*gamma(3/2)/pi^(3/2)*sqrt(b);
% mu=log(126)-s2/2;
% b = 1;
a = 1;
mu = 0;
rate = 2.6;

M    = 30; % number of realisations
Eps  = 2.^(-8:1:-5); % required accuracy

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

k = 4; % starting level

% MISMC TD
tic;
mi_LGC_com(@opre_miLGC_l,Eps,M,data,k,rate,a,b,mu,1);
toc;

% MISMC TP
tic;
mi_LGC_com(@opre_miLGC_l,Eps,M,data,k,rate,a,b,mu,2);
toc;