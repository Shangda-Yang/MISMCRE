% -------------------------------------------------------------------- %
% [sums, cost, u] = mlmc_l(l,N,theta,data)     
% single level routine
% inputs:  l     = level
%          N     = number of samples
%          sigma = std deviation of the error of y - G
%          data  = observations
% outputs: sums  = increments
%          cost  = computational cost
%          u     = samples generated from posterior distributions
% -------------------------------------------------------------------- %
function [sums, cost, u] = opre_mlmcmc_l(l,N,sigma,data,x_data)
y = data;

% generate a sample the from prior distribution
u_c = 2*rand(1) - 1;
u = zeros(1,N+1);
u(1) = u_c;
% solve pde at level l 
g_f = pde_solver(x_data, l, u(1)); 
% calculate the likelihoods
l_f_temp = exp(-0.5*(g_f - y)'*(g_f - y)/sigma/sigma); 

if l == 1
    g_c = 0;
    l_c_temp = 0;
else
    g_c = pde_solver(x_data, l-1, u(1));
    l_c_temp = exp(-0.5*(g_c - y)'*(g_c - y)/sigma/sigma);
end

% approximate coupling (maximal coupling)
G_k = zeros(1,N+1);
G_k(1) = max( l_f_temp, l_c_temp );

f_f = zeros(1,N+1);
f_c = zeros(1,N+1);

f_f(1) = l_f_temp;
f_c(1) = l_c_temp;

for k = 2:N+1
    u(k)   = u(k-1);
    G_k(k) = G_k(k-1);
    f_f(k) = f_f(k-1);
    f_c(k) = f_c(k-1);    
    % generate a sample from proposal distribution
    u_star = u(k-1) + 1.5*randn(1);
    if u_star <=1 && u_star >= -1
        % solve pde at level l and l-1
        g_f = pde_solver(x_data, l, u_star);
        
        % calculate the likelihoods
        l_f_temp = exp(-0.5*(g_f - y)'*(g_f - y)/sigma/sigma);
        
        if l == 1
            g_c =0;
            l_c_temp = 0;
        else
            g_c = pde_solver(x_data, l-1, u_star);
            l_c_temp = exp(-0.5*(g_c - y)'*(g_c - y)/sigma/sigma);
        end
        % approximate coupling (maximal coupling)
        G_star = max( l_f_temp, l_c_temp );
        % accept the sample with a probability of \alpha
        alpha_temp = G_star/G_k(k-1);
        alpha = min( 1, alpha_temp );
        uni = rand(1);
        if uni < alpha
            u(k) = u_star;
            G_k(k) = G_star;
            f_f(k) = l_f_temp;
            f_c(k) = l_c_temp;
        end
    end
end
% calculate the value of the s-n estimator
if l == 1
    sums = sum(u.^2)/N;    
    cost = 2^l;
else
    sums = sum(u.^2 .* f_f./G_k)/sum(f_f./G_k) - ...
        sum(u.^2 .* f_c./G_k)/sum(f_c./G_k);
    cost = 2*2^l;
end
% MCMC acceptance rate
% rate = sum(u(1:end-1) ~= u(2:end))/length(u);    
% % plot autocorrelation of samples
% figure(l)
% autocorr(u,30)
% shg
end

