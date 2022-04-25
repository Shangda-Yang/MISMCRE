% -------------------------------------------------------%
% function [sums, cost] = SMC(l,N,sigma,data,x_data)
% Single level SMC
% inputs:  l      = level
%          N      = number of paths
%          sigma  = std deviation of the error of y - G
%          data   = values of y
%          x_data = points corresponding to observations
%
% outputs: sums   = results of SMC
%          cost   = total computational costs
% -------------------------------------------------------%
function [sums, cost] = SMC(l,N,sigma,data,x_data)
y = data;
% set up the temporing step
lambda_diff = 0.5;
Lambda = (0:lambda_diff:1);

likelihood = zeros(1,N);

for j = 1:length(Lambda)
    lambda = Lambda(j); 
    if j == 1
        u = 2*rand(N,1)-1;
        H = ones(1,N);
    else
        % resampling step
        A = Multinomial_Resampling(W);
        u = u(A');
        % mutation step
        parfor k = 1:N
            % rate = 0;
            for i = 1:8
                u_star = u(k) + 0.5*randn(1);
                if u_star <= 1 && u_star >= -1
                    g_f = pde_solver(x_data, l, u_star);
                    
                    l_temp = exp(-0.5/sigma/sigma*(g_f - y)'*(g_f - y));
                    
                    alpha_temp = (l_temp/likelihood(k))^lambda;
                    
                    alpha = min( 1, alpha_temp );
                    
                    uni = rand(1);
                    
                    if uni < alpha
                        u(k)   = u_star;
                        likelihood(k) = l_temp;
                        H(k) = log( l_temp^lambda_diff );
                        % rate = rate + 1;
                    end
                end
            end
            % rate = rate/16;
        end
    end
    W = exp( H - min(H) )./sum(exp( H - min(H)) );
end
sums = sum(u.^2)/N;
cost = 2^l*N;
end
