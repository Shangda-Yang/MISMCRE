% ------------------------------------------------------- %
% [single_app,cost] = MCMC(l,N,theta,data,x_data)
% Single level MCMC
% inputs:  l      = level
%          N      = number of paths
%          sigma  = std deviation of the error of y - G
%          data   = observations
%          x_data = points corresponding to observations
%
% outputs: single_app = results of MCMC
%          cost       = total computational costs
% ------------------------------------------------------- %
function [single_app,cost] = MCMC(l,N,sigma,data,x_data)

y = data;
% generate a sample from prior distribution
u_c = 2*rand(1) - 1;

u = zeros(1,N+1);
u(1) = u_c;

g_f = pde_solver(x_data, l, u(1));
l_f_temp = exp(-0.5*(g_f - y)'*(g_f - y)/sigma/sigma); 

G_k = zeros(1,N+1);
G_k(1) = l_f_temp;

f_f = zeros(1,N+1);

f_f(1) = l_f_temp;

for k = 2:N+1
    u(k)   = u(k-1);
    G_k(k) = G_k(k-1);
    f_f(k) = f_f(k-1); 
    % generate a sample from proposal distribution
    u_star = u(k-1) + 0.8*randn(1);
    if u_star <=1 && u_star >= -1
    % solve pde
    g_f = pde_solver(x_data, l, u_star);
    
    % calculate the likelihood
    l_f_temp = exp(-0.5*(g_f - y)'*(g_f - y)/sigma/sigma);
    
    G_star = l_f_temp;
    
    % accept the sample with a probability of \alpha
    alpha_temp = G_star/G_k(k-1);
    
    alpha = min( 1, alpha_temp );
    
    uni = rand(1);
    
    if uni < alpha
        u(k) = u_star;
        G_k(k) = G_star;
        f_f(k) = l_f_temp;
    end
    end
    
end
single_app = sum(u.^2)/N;
cost = 2^l*N; 
end
