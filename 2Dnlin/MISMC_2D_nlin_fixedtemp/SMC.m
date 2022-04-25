% ------------------------------------------------------------------ %
% function mismc_complexity(mismc_l,Eps,sigma,data,x_data,k,setType)
%
% multi-index Monte Carlo complexity test routine
%
% input:  mlmc_l   = function for level l estimator
%         Eps      = vector of required accuracy
%         sigma    = std deviation of the error of y - G
%         data     = observations (values of y) vector
%         x_data   = corresponding observation points
%         k        = starting refinement level
%         setType  = type of index set
% .mat:   resn     = results of mimc with sn
%         rere     = results of mimc with re
%         costs    = total costs of mimc
%         mse_sn   = mse of self-normalised increments estimator
%         mse_re   = mse of ratio estimator
%         rLf      = actural finest level of refinement
% ------------------------------------------------------------------ %
function SMC(SMCL,Eps,sigma,data,x_data)
% number of realisations
rel = 50;
% results
re = zeros(rel,length(Eps));
% costs
costs = zeros(1,length(Eps));
% loading emperical results first
load('./workingdata/mlsmc_emperical.mat','se1','cost')
load('./workingdata/solutions.mat','solutions')

for i = 1:length(Eps)
    eps = Eps(i);
    L = ceil(log2(1/eps/3)/2);
    N = ceil(se1(1)/eps/eps);
    costl = 0;
    parfor j = 1:rel    
        [sums, cst] = SMCL(L,N,sigma,data,x_data);
        re(j,i) = sums;
        costl = costl + cst;
    end
    costs(i) = costl/rel;
end

solution1 = solutions(1);
solution2 = solutions(2);

mse1 = sum((re - solution1).^2)./rel;
mse2 = sum((re - solution2).^2)./rel;

save('./workingdata/smc_complexity.mat','re',...
    'mse1','mse2','costs','L')
end