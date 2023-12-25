% --------------------------------------------------------- %
% function mlmcmc_test(mlmcmc_l, N, L, Eps, fp, sigma, data, x_data)
%
% multilevel Markov Chain Monte Carlo test routine
%
% mlmcmc_l = function for level l estimator
% N        = number of samples for convergence tests
% L        = number of levels for convergence tests
% Eps      = desired accuracy array for MLMC calcs
% fp       = file handle for printing to file
% sigma    = std deviation of the error of y - G
% data     = observations
% x_data   = points corresponding to observations
% --------------------------------------------------------- %
function mlmcmc_test(mlmcmc_l, N, L, Eps, fp, sigma, data, x_data)

%
% first, convergence tests
%

PRINTF2(fp,'\n');
PRINTF2(fp,'**********************************************************\n');
PRINTF2(fp,'*** MATLAB mlmcmc_test on %s         ***\n',datestr(now));
PRINTF2(fp,'**********************************************************\n');
PRINTF2(fp,'\n');
PRINTF2(fp,'*** Convergence tests                                  ***\n');
PRINTF2(fp,'*** using N =%7d samples                           ***\n',N);
PRINTF2(fp,'**********************************************************\n');
PRINTF2(fp,'\n l   sn_inc(phi)    inc_var   cst');
PRINTF2(fp,'\n--------------------');
PRINTF2(fp,'----------------------\n');

del1 = [];
var1 = [];
cost = [];

for l = 1:L  

    sums = [];
    cst  = 0;
    
    M = 100;
    
    parfor j = 1:M   
        [sums_j, cst_j, ~] = mlmcmc_l(l,N,sigma,data,x_data);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;
    end
    
        cost = [cost cst];
        del1 = [del1 abs(sum(sums)/M)];
        var1 = [var1 N*sum( (sums - sum(sums)/M).^2 )/M ];
    
%     PRINTF2(fp,'%2d  %11.4e  %.3e  %.2e\n', ...
%         l,del1(l),var1(l),cst);
    PRINTF2(fp,'%2d  %11.4e   %.3e  %.2e \n', ...
        l,del1(l),var1(l),cst);

end

%
% use linear regression to estimate alpha, beta and gamma, and c1, c2 and c3
%

L1 = 2;
L2 = L;
pa    = polyfit(L1:L2,log2(abs(del1(L1:L2))),1);  alpha = -pa(1); c1 = 2^pa(2);
pb    = polyfit(L1:L2,log2(abs(var1(L1:L2))),1);  beta  = -pb(1); c2 = 2^pb(2);
pg    = polyfit(L1:L2,log2(abs(cost(L1:L2))),1);  gamma =  pg(1); c3 = 2^pg(2);

PRINTF2(fp,'\n********************************************************\n');
PRINTF2(fp,'*** Linear regression estimates of MLMCMC parameters ***\n');
PRINTF2(fp,'********************************************************\n');
PRINTF2(fp,'\n alpha = %f  (exponent for MLMC weak convergence)\n',alpha);
PRINTF2(fp,' beta  = %f  (exponent for MLMC variance) \n',beta);
PRINTF2(fp,' gamma = %f  (exponent for MLMC cost) \n',gamma);
PRINTF2(fp,' c1  = %f  (constant for MLMC weak convergence)\n',c1);
PRINTF2(fp,' c2  = %f  (constant for MLMC variance) \n',c2);
PRINTF2(fp,' c3  = %f  (constant for MLMC cost) \n',c3);


%
% second, mlmc complexity tests
%

PRINTF2(fp,'\n');
PRINTF2(fp,'******************************* \n');
PRINTF2(fp,'*** MLMCMC complexity tests *** \n');
PRINTF2(fp,'******************************* \n\n');
PRINTF2(fp,'   eps    L    mlmc_cost   slmc_cost   mse_ml   mse_sl   N_l\n');
PRINTF2(fp,'-----------------------------------------------------------------------\n');


alpha = 2;
beta  = 4;

% number of realisations for MSE
M = 100;

% analytic solution of the quantity of interest (u^2)
u2 = analytic_u2(sigma,data,x_data);

for i = 1:length(Eps)
    eps = Eps(i);
    
    Cl   = [];
    slCl = [];
    mse_ml = 0;
    mse_sl = 0;
    
%%%%% number of levels required for eps by constants
    deno = (2^alpha - 1)*eps;
    L = ceil( log2(sqrt(2)*c1/deno)/alpha );    
    
%%%%% number of samples required for eps by constants
%     L_set = (1:1:L);
%     K_l = sqrt(var1(1)*2^(-1))+ ...
%         sum( sqrt(c2*c3*2.^(-L_set(2:end)*(beta - gamma))) );
%         
%     N_l = ceil(  2*eps^(-2)*K_l.*sqrt(c2*2.^-(L_set.*(beta+gamma)))/c3);
%     
%     N_l(1) = ceil( 2*eps^(-2)*K_l*sqrt(var1(1)*2) );

    N_l = ceil( 2*eps^(-2)*sqrt(sum(var1(1:L).*cost(1:L))).* ...
        sqrt(var1(1:L)./cost(1:L)) );

    % calculate MSE
    parfor j = 1:M
        [P,cl] = mlmcmc(mlmcmc_l,eps,L,N_l,sigma,data,x_data);
        mse_ml = mse_ml + (P - u2)^2/M;
        Cl = [Cl; cl./M];
    end
    
    deno = (2^alpha - 1)*eps;
    sL = max(1,ceil( log2(sqrt(2)/deno)/alpha )); 
    slN = ceil(2*var1(1)/eps/eps);
    parfor j = 1:M
        [slP, slcl] = MCMC(sL,slN,sigma,data,x_data);
        mse_sl = mse_sl + (slP - u2)^2/M;
        slCl = [slCl; slcl/M];
    end
    
    mlmc_cost = sum(N_l.*sum(Cl));
    slmc_cost = sum(slCl);
    PRINTF2(fp,'%.3e %d %.3e  %.3e  %.3e  %.3e', ...
        eps,L, mlmc_cost,slmc_cost, mse_ml,mse_sl);
    PRINTF2(fp,'%9d',N_l);
    PRINTF2(fp,'\n');

end

PRINTF2(fp,'\n');

end

%
% function to print to both a file and stdout 
%

function PRINTF2(fp,varargin)
fprintf(fp,varargin{:});
fprintf( 1,varargin{:});
end
