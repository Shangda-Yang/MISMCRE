% ------------------------------------------------------------------ %
% function mismc_emperical(mimc_l,N, Lx, Ly,sigma,data)
%
% multi-index Monte Carlo mixed regularity test routine
% and emperical mean and variance calculation
%
% input:   mlsmc_l  = function for level l estimator
%          N        = number of samples for convergence tests
%          M        = number of realisations
%          Lmax     = max number of levels for convergence tests in x and y
%          data     = observation points
%          k        = starting level
%          rate     = \beta^{\prime}
%          a        = \theta_2
%          b        = \theta_3
%          mu       = \theta_1
% outputs: we{1}    = difference with self-normalised increments estimator
%          we{2}    = difference unnormalised integral
%          we{3}    = difference normalising constant
%          se{1}    = variance of sn increments estimtor
%          se{2}    = unnormalised integral of squared target
%          se{3}    = unnormalised integral of 1
%          cost     = cost at each level
% ------------------------------------------------------------------ %

function [we,se,cost] = mi_LGC_emp(mismc_l,N,M,Lmax,data,k,rate,a,b,mu)
we1  = zeros(Lmax-k+1,Lmax-k+1);
we2  = zeros(Lmax-k+1,Lmax-k+1);
we3  = zeros(Lmax-k+1,Lmax-k+1);
se1  = zeros(Lmax-k+1,Lmax-k+1);
se2  = zeros(Lmax-k+1,Lmax-k+1);
se3  = zeros(Lmax-k+1,Lmax-k+1);
se4  = zeros(Lmax-k+1,Lmax-k+1);
se5  = zeros(Lmax-k+1,Lmax-k+1);
se6  = zeros(Lmax-k+1,Lmax-k+1);
cost = zeros(Lmax-k+1,Lmax-k+1);

for lx = k:Lmax
    for ly = k:Lmax
    sums = [];
    cst  = 0;
    
    fprintf(1,'lx = %.1d, ly = %.1d \n',lx,ly)
    
    parfor j = 1:M
        [sums_j, cst_j] = mismc_l(lx,ly,N,data,rate,a,b,mu,k);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;      
    end
    
    idx = lx - k + 1;
    idy = ly - k + 1;
    
    cost(idx,idy)  = cst;
    we1(idx,idy)  = abs( sum(sums(:,1))/M );
    we2(idx,idy)  = abs( sum(sums(:,2))/M );
    we3(idx,idy)  = abs( sum(sums(:,3))/M );
    
    se1(idx,idy)  = N*sum( (sums(:,1) - sum(sums(:,1))/M).^2 )/M;
    se2(idx,idy)  = N*sum( (sums(:,2) - sum(sums(:,2))/M).^2 )/M;
    se3(idx,idy)  = N*sum( (sums(:,3) - sum(sums(:,3))/M).^2 )/M;
    
    se4(idx,idy)  = sum(sums(:,1).^2)/M;
    se5(idx,idy)  = sum(sums(:,2).^2)/M;
    se6(idx,idy)  = sum(sums(:,3).^2)/M;

    end
end
we = {we1,we2,we3};
se = {se1,se2,se3,se4,se5,se6};
end

