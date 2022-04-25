% ----------------------------------------------------------- %
% function [sums, cost] = opre_mlsmc_l(l,N,sigma,data,x_data)
% single level routine
% inputs:  l      = level
%          N      = number of samples
%          sigma  = std deviation of the error of y - G
%          data   = values of y
%          x_data = points corresponding to observations
% outputs: sums(1) = SN increments estimator at level l
%          sums(2) = unnormalised increments phi at level l
%          sums(3) = unnormalised increments 1 at level l
%          cost    = computational cost
% ----------------------------------------------------------- %
function [sums, cost] = opre_mlsmc_l(l,N,sigma,data,x_data)
y = data;
sums(1:3) = 0;

% set up the temporing step
lambda_diff = 0.5;
Lambda = (lambda_diff:lambda_diff:1);

% initialisation
u = 2*rand(N,1)-1;
Z = 1;

% multi-level increments
H   = zeros(N,1);
G_k = zeros(N,1);
f_f = zeros(N,1);
f_c = zeros(N,1);

parfor i = 1:N
    g_f = pde_solver(x_data, l, u(i));
    l_f_temp = exp(-0.5*(g_f - y)'*(g_f - y)/sigma/sigma);
    if l == 1
        g_c = 0;
        l_c_temp = 0;
    else
        g_c = pde_solver(x_data, l-1, u(i));
        l_c_temp = exp(-0.5*(g_c - y)'*(g_c - y)/sigma/sigma);
    end
    G_k(i) = max( l_f_temp, l_c_temp );
    f_f(i) = l_f_temp;
    f_c(i) = l_c_temp;
    H(i) = log( max( l_f_temp, l_c_temp )^lambda_diff );
end


for j = 1:length(Lambda)
    
    lambda = Lambda(j);
    % normalising constant
    Z = Z*sum(exp(H))/N;
    if j ~= 1
        W = exp( H - min(H) )./sum(exp( H - min(H)) );
        % fprintf(1,'ESS: %d \n',1/sum(W.^2));
        % resampling step
        A = Multinomial_Resampling(W);
        u = u(A');
    end
    parfor k = 1:N
%                        rate = 0;
        for i = 1:9
            u_star = u(k) + 0.5*randn(1);
            if u_star <= 1 && u_star >= -1 
                g_f = pde_solver(x_data, l, u_star);
                l_f_temp = exp(-0.5*(g_f - y)'*(g_f - y)/sigma/sigma);
                if l == 1
                    g_c =0;
                    l_c_temp = 0;
                else
                    g_c = pde_solver(x_data, l-1, u_star);
                    l_c_temp = exp(-0.5*(g_c - y)'*(g_c - y)/sigma/sigma);
                end
                G_star = max( l_f_temp, l_c_temp );
                alpha_temp = (G_star/G_k(k))^lambda;
                alpha = min( 1, alpha_temp );
                uni = rand(1);
                if uni < alpha
                    u(k) = u_star;
                    G_k(k) = G_star;
                    
                    f_f(k) = l_f_temp;
                    f_c(k) = l_c_temp;
                    H(k) = log( G_star^lambda_diff );
%                       rate = rate + 1;
                end
            end
        end
%           u_rate = rate/8;
%          fprintf(1,'%.5f \n',u_rate)
    end
end

if l == 1
    sums(1) = sum(u.^2)/N;
    sums(2) = Z*sum(u.^2)/N;
    sums(3) = Z;
    cost = 2^l;
else
    sums(1) = sum(u.^2 .* f_f./G_k)/sum(f_f./G_k) - ...
        sum(u.^2 .* f_c./G_k)/sum(f_c./G_k);
    sums(2) = Z*( sum(u.^2 .* f_f./G_k)- sum(u.^2 .* f_c./G_k) )/N;
    sums(3) = Z*( sum(f_f./G_k) - sum(f_c./G_k) )/N;
    cost = 2*2^l;
end

end
