% ------------------------------------------------------------------ %
% function mlsmc_complexity(mismc_l,Eps,sigma,data,x_data,k,setType)
%
% Multilevel Monte Carlo complexity test routine
%
% input:  mlsmc_l   = function for level l estimator
%         Eps      = vector of required accuracy
%         sigma    = std deviation of the error of y - G
%         data     = observations (values of y) vector
%         x_data   = corresponding observation points
%         k        = starting refinement level
% .mat:   resn     = results of mimc with sn
%         rere     = results of mimc with re
%         costs    = total costs of mimc
%         mse_sn   = mse of self-normalised increments estimator
%         mse_re   = mse of ratio estimator
% ------------------------------------------------------------------ %

function mlsmc_complexity(mlsmc_l,Eps,sigma,data,x_data,k)
[filepart,~,~] = fileparts(pwd);
% number of realisations
rel = 50;
% results
resn = zeros(rel,length(Eps));
rere = zeros(rel,length(Eps));
% costs
costs = zeros(1,length(Eps));
% loading emperical results first
loadpath = fullfile(filepart, 'Results','mlsmc','emperical','mlsmc_emperical.mat');
load(loadpath,'se1','cost')
loadpath = fullfile(filepart, 'Results','solutions.mat');
load(loadpath,'solutions')

for i = 1:length(Eps)
    eps = Eps(i);
    L = ceil(log2(sqrt(2)/eps/3)/2)-k;
    vc1 = sum(sum(sqrt(se1(1:L+1).*cost(1:L+1)))).*...
        sqrt(se1(1:L+1)./cost(1:L+1));
    M = 2*ceil(vc1./eps/eps);
    costl = 0;
    parfor j = 1:rel    
        sumlsn  = 0;
        sumlret = 0;
        sumlreb = 0;
        % total degree index set
        for l = 0:L
            N = M(l+1);
            [sums, cst] = mlsmc_l(l,l,N,sigma,data,x_data,k);
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

solution1 = solutions(1);
solution2 = solutions(2);

mse_sn1 = sum((resn - solution1).^2)./rel;
mse_re1 = sum((rere - solution1).^2)./rel;
mse_sn2 = sum((resn - solution2).^2)./rel;
mse_re2 = sum((rere - solution2).^2)./rel;

rL = L + k;
savepath = fullfile(filepart, 'Results','mlsmc','complexity','mlsmc_complexity.mat');
save(savepath,'resn','rere',...
    'mse_sn1','mse_re1','mse_sn2','mse_re2','costs','rL','k')

end
