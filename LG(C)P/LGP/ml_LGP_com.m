% ---------------------------------------------------------------- %
% function ml_LGC_com(mlsmc_l,Eps,rel,data,k,rate,a,b,mu)
% inputs:   mlsmc_l = mlsmc single level routine
%           Eps     = required accuracy
%           rel     = number of realisations
%           k       = starting level
%           rate    = \beta^{\prime}
%           a       = \theta_2
%           b       = \theta_3
%           mu      = \theta_1
% outputs:  mat file with data
% ---------------------------------------------------------------- %

function ml_LGP_com(mlsmc_l,Eps,rel,data,k,rate,a,b,mu)
% results
resn = zeros(rel,length(Eps));
rere = zeros(rel,length(Eps));
% costs
costs = zeros(1,length(Eps));
% loading emperical results first
% load('./workingdata/mlsmc_emperical.mat','se1','cost')
% load('./workingdata/solution.mat','solution')

for i = 1:length(Eps)
    eps = Eps(i)
    alpha = 0.8;
    beta  = 1.6;
    L = max(ceil(log2(1/eps/(2^alpha-1))/alpha),k+1);
    L_set = (k:1:L);
    var   = 2.^(-beta*L_set);
    var(1) = 1e-3;
    cost  = 2*L_set.*2.^(2*L_set);
    K_l = sum(sum( sqrt(var.*cost) ));
    M = max(ceil( eps^(-2)*K_l.*sqrt(var./cost) ), 2);
    costl = 0;
    parfor j = 1:rel
        sumlsn  = 0;
        sumlret = 0;
        sumlreb = 0;
        % total degree index set
        for l = k:L
            N = M(l-k+1);
            [sums, cst] = mlsmc_l(l,l,N,data,rate,a,b,mu,k);
            sumlsn  = sumlsn  + sums(1);
            sumlret = sumlret + sums(2);
            sumlreb = sumlreb + sums(3);
            costl   = costl + N*cst;
        end
        resn(j,i) = sumlsn;
        rere(j,i) = sumlret/sumlreb;
    end
    costs(i) = costl/rel;
end

% mse_sn = sum((resn - solution).^2)./rel;
% mse_re = sum((rere - solution).^2)./rel;

sn = resn;
re = rere;

% save('./workingdata/mismc_complexity_td.mat','resn','rere',...
%     'mse_sn','mse_re','costs','L','k')
save('./workingdata/mlcomplexity.mat','sn','re','costs','L','k')

end



