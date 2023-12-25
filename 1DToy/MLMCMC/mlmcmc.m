% ------------------------------------------------------------ %
% function [P, Cl] = mlmcmc(mlmcmc_l,eps,L,N_l,sigma,data,x_data)
% multi-level Monte Carlo estimation
% inputs:  mlmc_l = function for level l estimator
%          eps    = desired accuracy (rms error)           > 0
%          L      = number of levels for convergence tests > 0
%          N_l    = number of samples at different levels  > 0
%          sigma  = std deviation of the error of y - G
%          data   = observations
%          x_data = points corresponding to observations
% outputs: P      = E[quantity of interest]
%          cost   = total cost
% ------------------------------------------------------------ %

function [P, Cl] = mlmcmc(mlmcmc_l,eps,L,N_l,sigma,data,x_data)

%
% check input parameters
%

if (L<=0)
    error('error: needs L > 0 \n');
end

if (N_l<=0)
    error('error: needs N_l>0 \n');
end

if (eps<=0)
    error('error: needs eps>0 \n');
end

suml(1,1:L) = 0;
costl(1:L)    = 0;

% update sample sums

for l=1:L
    Nl = N_l(l);
    [sums, cost, ~] = mlmcmc_l(l,Nl,sigma,data,x_data);
    suml(1,l) = sums;
    costl(l)  = cost;
end

% evaluate multilevel estimator and cost
P = sum(suml(1,:));
Cl = costl;

end
