% ------------------------------------------------------------------ %
% function mismc_emperical(mimc_l,N, Lx, Ly,sigma,data,x_data)
%
% multi-index Monte Carlo assumption test routine
% and emperical mean and variance calculation
%
% input:  mlmc_l   = function for level l estimator
%         N        = number of samples for convergence tests
%         Lx       = number of levels for convergence tests in x
%         Ly       = number of levels for convergence tests in y
%         sigma    = std deviation of the error of y - G
%         data     = observations (values of y)
%         x_data   = corresponding observation points
% .mat:   we1      = difference with self-normalised increments estimator
%         we2      = difference unnormalised integral
%         we3      = difference normalising constant
%         se1      = variance of sn increments estimtor
%         se2      = unnormalised integral of squared target
%         se3      = unnormalised integral of 1
%         cost     = cost at each level
% ------------------------------------------------------------------ %
function mismc_emperical(mismc_l,N, Lx, Ly,sigma,data,x_data,K)

we1  = zeros(Lx+1,Ly+1);
we2  = zeros(Lx+1,Ly+1);
we3  = zeros(Lx+1,Ly+1);

se1  = zeros(Lx+1,Ly+1);
se2  = zeros(Lx+1,Ly+1);
se3  = zeros(Lx+1,Ly+1);

cost = zeros(Lx+1,Ly+1);

for lx = 0:Lx
    for ly = 0:Ly
    
    sums = [];
    cst  = 0;
    
    M = 20;
    
    parfor j = 1:M  
        [sums_j, cst_j] = mismc_l(lx,ly,N,sigma,data,x_data,K);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;   
    end
    
    idx = lx + 1;
    idy = ly + 1;
    
    cost(idx,idy) = cst;
    we1(idx,idy)  = abs( sum(sums(:,1))/M );
    we2(idx,idy)  = abs( sum(sums(:,2))/M );
    we3(idx,idy)  = abs( sum(sums(:,3))/M );
    
    se1(idx,idy)  = N*sum( (sums(:,1) - sum(sums(:,1))/M).^2 )/M;
    se2(idx,idy)  = N*sum( (sums(:,2) - sum(sums(:,2))/M).^2 )/M;
    se3(idx,idy)  = N*sum( (sums(:,3) - sum(sums(:,3))/M).^2 )/M;
    end
end

%
% use linear regression to estimate alpha, beta and gamma, and c1, c2 and c3
%
% 
% pa  = polyfit(1+k:Lx+k,log2(abs(we1(2:end,Ly+1))),1);  alphax = -pa(1); c1x = 2^pa(2);
% pb  = polyfit(1+k:Lx+k,log2(abs(se1(2:end,Ly+1))),1);  betax  = -pb(1); c2x = 2^pb(2);
% pg  = polyfit(1+k:Lx+k,log2(abs(cost(2:end,Ly+1))),1); gammax =  pg(1); c3x = 2^pg(2);
% 
% pa  = polyfit(1+k:Ly+k,log2(abs(we1(Lx+1,2:end))),1);   alphay = -pa(1); c1y = 2^pa(2);
% pb  = polyfit(1+k:Ly+k,log2(abs(se1(Lx+1,2:end))),1);   betay  = -pb(1); c2y = 2^pb(2);
% pg  = polyfit(1+k:Ly+k,log2(abs(cost(Lx+1,2:end))),1);  gammay =  pg(1); c3y = 2^pg(2);
% 
% well = diag(we1);
% sell = diag(se1);
% costll = diag(cost);
% 
% ll = ( [0:Lx]+k ).*2;
% 
% pa  = polyfit(ll(2:end),log2(abs(well(2:end))),1);    alpha = -pa(1); c1 = 2^pa(2);
% pb  = polyfit(ll(2:end),log2(abs(sell(2:end))),1);    beta  = -pb(1); c2 = 2^pb(2);
% pg  = polyfit(ll(2:end),log2(abs(costll(2:end))),1);  gamma =  pg(1); c3 = 2^pg(2);
% 
% fprintf(1,'alphax = %.4f',alphax)
% fprintf(1,'betax  = %.4f',betax)
% fprintf(1,'gammax = %.4f',gammax)
% fprintf(1,'alphay = %.4f',alphay)
% fprintf(1,'betay  = %.4f',betay)
% fprintf(1,'gammay = %.4f',gammay)
% fprintf(1,'alpha  = %.4f',alpha)
% fprintf(1,'beta   = %.4f',beta)
% fprintf(1,'gamma  = %.4f',gamma)

savepath = fullfile(filepart, 'Results','mismc','emperical','mismc_emperical.mat');
save(savepath,'Lx','Ly','K',...
    'cost','we1','we2','we3','se1','se2','se3')...
%     'alphax','betax','gammax','alphay','betay','gammay',...
%     'alpha','beta','gamma',...
%     'c1x','c1y','c1','c2x','c2y','c2','c3x','c3y','c3')
end

