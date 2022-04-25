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
function mismc_complexity(mismc_l,Eps,sigma,data,x_data,K,setType)
% number of realisations
rel = 50;
% results
resn = zeros(rel,length(Eps));
rere = zeros(rel,length(Eps));
% costs
costs = zeros(1,length(Eps));
% loading emperical results first
loadpath = fullfile(filepart, 'Results','mismc','emperical','mismc_emperical.mat');
load(loadpath,'se1','cost')
loadpath = fullfile(filepart, 'Results','solutions.mat');
load(loadpath,'solutions')

if setType == 1
for i = 1:length(Eps)
    eps = Eps(i);
    L = ceil(log2(2/eps)/2)-K;
    vc1 = sum(sum(sqrt(se1(1:L+1,1:L+1).*cost(1:L+1,1:L+1)))).*...
        sqrt(se1(1:L+1,1:L+1)./cost(1:L+1,1:L+1));
    M = ceil(vc1./eps/eps);
    costl = 0;
    parfor j = 1:rel    
        sumlsn  = 0;
        sumlret = 0;
        sumlreb = 0;
        % total degree index set
        for lx = 0:L
            for ly = 0:L
                N = M(lx+1,ly+1);
                [sums, cst] = mismc_l(lx,ly,N,sigma,data,x_data,K);
                sumlsn  = sumlsn  + sums(1);
                sumlret = sumlret + sums(2);
                sumlreb = sumlreb + sums(3);
                costl   = costl + N*cst;
            end
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

rL = L + K;
savepath = fullfile(filepart, 'Results','mismc','complexity','mismc_complexity_tp.mat');
save(savepath,'resn','rere',...
    'mse_sn1','mse_re1','mse_sn2','mse_re2','costs','rL','K')  
elseif setType == 2
    
for i = 1:length(Eps)
    eps = Eps(i);
    L = ceil( 0.25/log(2)*(log((log(1/eps))^2/eps)) )-K;
    se = fliplr(triu(fliplr(se1(1:L+1,1:L+1)))); 
    vc1 = sum(sum(sqrt(se.*cost(1:L+1,1:L+1)))).*...
        sqrt(se(1:L+1,1:L+1)./cost(1:L+1,1:L+1));
    M = ceil( vc1./eps/eps );
    costl = 0;
    parfor j = 1:rel    
        sumlsn  = 0;
        sumlret = 0;
        sumlreb = 0;
        % total degree index set
        for lx = 0:L
            for ly = 0:L-lx
                N = M(lx+1,ly+1);
                [sums, cst] = mismc_l(lx,ly,N,sigma,data,x_data,K);
                sumlsn  = sumlsn  + sums(1);
                sumlret = sumlret + sums(2);
                sumlreb = sumlreb + sums(3);
                costl   = costl + N*cst;
            end
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

rL = L + K;
savepath = fullfile(filepart, 'Results','mismc','complexity','mismc_complexity_td.mat');
save(savepath,'resn','rere',...
    'mse_sn1','mse_re1','mse_sn2','mse_re2','costs','rL','K')    
else
    fprintf(1,'Wrong setType')
end
end

