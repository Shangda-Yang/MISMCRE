close all; clear all;

randn('seed',20)
rand('seed',20)

s2=1.91;
b=(33/pi)^2;
% a=s2*gamma(3/2)/pi^(3/2)*sqrt(b);
% mu=log(126)-s2/2;
% b = 1;
a = 1;
mu = 0;

M    = 100;
Eps  = 2.^(-8:1:-3);

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

rate = 2.6;

k = 5;
tic;
ml_LGC_com(@opre_mlLGC_l,Eps,M,data,k,rate,a,b,mu);
toc;