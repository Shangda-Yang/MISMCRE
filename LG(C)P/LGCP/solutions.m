close all; clear all;
% for reference solutions
randn('seed',20)
rand('seed',20)

s2=1.91;
b=(33/pi)^2;
% a=s2*gamma(3/2)/pi^(3/2)*sqrt(b);
% mu=log(126)-s2/2;
% b = 1;
a = 1;
mu = 0;

eps = 2^-9;

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

rate = 2.6;

k = 4;
tic;
[solution] = getSol(@opre_mlLGC_l,eps,data,k,rate,a,b,mu);
toc;

save('./workingdata/solutions.mat','solution')

function [solution] = getSol(mlsmc_l,eps,data,k,rate,a,b,mu)
% load('./workingdata/mlsmc_emperical.mat','se1','cost')

alpha = 0.8;
beta  = 1.6;
L = max(ceil(log2(1/eps/(2^alpha-1))/alpha),k+1);
L_set = (k:1:L);
var   = 2.^(-beta*L_set);
cost  = 2.^(2*L_set);
K_l = sum(sum( sqrt(var.*cost) ));
M = ceil(  eps^(-2)*K_l.*sqrt(var./cost) );

sumlsn  = 0;
sumlret = 0;
sumlreb = 0;
costl   = 0;

% total degree index set
for l = k:L
    N = M(l-k+1);
    [sums, cst] = mlsmc_l(l,l,N,data,rate,a,b,mu,k);
    sumlsn  = sumlsn  + sums(1);
    sumlret = sumlret + sums(2);
    sumlreb = sumlreb + sums(3);
    costl   = costl + N*cst;
end

solution(1) = sumlsn;
solution(2) = sumlret/sumlreb;

end