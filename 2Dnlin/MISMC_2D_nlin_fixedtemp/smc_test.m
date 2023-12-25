close all; clear all;

randn('seed',20)
rand('seed',20)

sigma  = 0.5;
 
Eps = [ 0.001 0.0025 0.005 0.01 0.025 ];
% observation points
x_data = [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75]; 
load('./workingdata/observations.mat','data')

tic;
SMC(@opre_SMC,Eps,sigma,data,x_data);
toc;

% ------------------------------------------------------------------ %
% [sums, cost] = opre_SMC(l,N,sigma,data,x_data)   
% SMC routine
%
% inputs:  l      = level of refinment
%          N      = number of samples
%          sigma  = std deviation of the error of y - G
%          data   = observations
%          x_data = corresponding points of observations
%
% output:  sums    = quantity of interests
%          cost    = computational cost
% ------------------------------------------------------------------ %
function [sums, cost] = opre_SMC(l,N,sigma,data,x_data)
y = data;

% set up the temporing step
lambda_diff = 0.5;
Lambda = (lambda_diff:lambda_diff:1);

% initialisation
u_single = 2*rand(N,2)-1; 
Z_single = 1;

% single level approximation
H_single = zeros(1,N);
f_single = zeros(1,N);
parfor i = 1:N
    [g_single, ~] = pde_solver_2D(l,l,x_data,u_single(i,:));
    l_f_temp = exp(-0.5/sigma/sigma*(g_single - y)'*(g_single - y));
    f_single(i) = l_f_temp;
    H_single(i) = log(l_f_temp^lambda_diff);
end

for j = 1:length(Lambda)
    
    lambda = Lambda(j);
    
    % normalising constant
    Z_single = Z_single*sum(exp(H_single))/N;
    
    W = exp(H_single - min(H_single))./ ...
        sum(exp(H_single - min(H_single)));
    %
    %       fprintf(1,'ESS: %d \n',1/sum(W.^2));
    
    % resampling step
    A = Multinomial_Resampling(W);
    
    u_single = u_single(A',:);
    
    % mutation step
    
    parfor k = 1:N
%          rate = 0;
        for i = 1:8
            u_star = u_single(k,:) + 1.5*randn(1,2);
            
            [g_single, ~] = pde_solver_2D(l,l,x_data,u_star);
            
            l_f_temp = exp(-0.5/sigma/sigma*(g_single - y)'*(g_single - y));
            
            if u_star(1) >= 1 || u_star(1) <= -1 || u_star(2) >= 1 || u_star(2) <= -1
                alpha_temp = 0;
            else
                alpha_temp = (l_f_temp/f_single(k))^lambda;
            end
            
            alpha = min( 1, alpha_temp );
            
            uni = rand(1);
            
            if uni < alpha
                u_single(k,:) = u_star;
                f_single(k) = l_f_temp;
                H_single(k) = log(l_f_temp^lambda_diff);
%                  rate = rate + 1;
            end
        end
%         u_rate = rate/1000;
%         fprintf(1,'acc_rate: %.4f \n',u_rate);
%         figure(1)
%         autocorr(u_single(:,2),30)
%         shg
    end
end
L2norm2 = u_single(:,1).^2 + u_single(:,2).^2;
sums = sum(L2norm2)/N;

cost = 2^(2*l);
end
