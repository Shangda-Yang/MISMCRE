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

function [we,se,cost] = ml_LGC_emp(mlsmc_l,N,M,Lmax,data,k,rate,a,b,mu)
we   = zeros(3,Lmax-k+1);
se   = zeros(3,Lmax-k+1);
cost = zeros(1,Lmax-k+1);

for l = k:Lmax
    fprintf(1,'Level %.1d \n',l)
    sums = [];
    cst  = 0;
    
    for j = 1:M
        [sums_j, cst_j] = mlsmc_l(l,l,N,data,rate,a,b,mu,k);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;      
    end
  
    id = l - k + 1;
    
    cost(id)  = cst;
    we(1,id)  = abs( sum(sums(:,1))/M );
    we(2,id)  = abs( sum(sums(:,2))/M );
    we(3,id)  = abs( sum(sums(:,3))/M );
    
    se(1,id)  = N*sum( (sums(:,1) - sum(sums(:,1))/M).^2 )/M;
    se(2,id)  = N*sum( (sums(:,2) - sum(sums(:,2))/M).^2 )/M;
    se(3,id)  = N*sum( (sums(:,3) - sum(sums(:,3))/M).^2 )/M;

end

end

