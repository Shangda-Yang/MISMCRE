close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver2Dnlin'));

N = 2000;
M = 20;
sigma  = 0.5;
% observation points
x_data = [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75];
loadpath = fullfile(filepart, 'Results','observations.mat');
load(loadpath,'data')

k = 2;

L = 4;

we1 = zeros(1,L+1);
we2 = zeros(1,L+1);
we3 = zeros(1,L+1);

se1 = zeros(1,L+1);
se2 = zeros(1,L+1);
se3 = zeros(1,L+1);

cost = zeros(1,L+1);

tic;
for l = 0:L
    
    sums = [];
    cst  = 0;
    parfor j = 1:M
        [sums_j, cst_j] = opre_MLSMC_l(l,l,N,sigma,data,x_data,k);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;  
    end
    id = l + 1;
    cost(1,id) = cst;
    
    we1(1,id)  = abs( sum(sums(:,1))/M );
    we2(1,id)  = abs( sum(sums(:,2))/M );
    we3(1,id)  = abs( sum(sums(:,3))/M );
    
    se1(1,id)  = N*sum( (sums(:,1) - sum(sums(:,1))/M).^2 )/M;
    se2(1,id)  = N*sum( (sums(:,2) - sum(sums(:,2))/M).^2 )/M;
    se3(1,id)  = N*sum( (sums(:,3) - sum(sums(:,3))/M).^2 )/M;
end
toc;
savepath = fullfile(filepart, 'Results','mlsmc','emperical','mlsmc_emperical.mat');
save(savepath,'k','L',...
    'we1','we2','we3','se1','se2','se3','cost')
rmpath(fullfile(filepart,'solver2Dnlin'));
