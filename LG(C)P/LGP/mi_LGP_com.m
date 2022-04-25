% ---------------------------------------------------------------- %
% function mi_LGP_com(mismc_l,Eps,rel,data,k,rate,a,b,mu,setType)
% inputs:   mismc_l = mismc single level routine
%           Eps     = required accuracy
%           rel     = number of realisations
%           k       = starting level
%           rate    = \beta^{\prime}
%           a       = \theta_2
%           b       = \theta_3
%           mu      = \theta_1
%           setType = 1: TD;2: TP
% outputs:  mat file with data
% ---------------------------------------------------------------- %

function mi_LGP_com(mismc_l,Eps,rel,data,k,rate,a,b,mu,setType)
% results
resn = zeros(rel,length(Eps));
rere = zeros(rel,length(Eps));
% costs
costs = zeros(1,length(Eps));
% loading emperical results first
% load('./workingdata/mismc_emperical.mat','se1','cost')
% load('./workingdata/solution.mat','solution')
switch setType
    case 1
        for i = 1:length(Eps)
            eps = Eps(i);
            alpha = 0.8;
            beta  = 1.6;
            
%             temp = 0.5/log(2)/alpha;
%             logEps = -log(eps);
%             L = ceil(temp*(logEps + log(temp*logEps) + ...
%                 4*exp(log(2)*alpha*sqrt(2))*(temp^2+temp) ));
            L = max(ceil(log((log(1/eps))^2/eps)/log(2)/2/alpha), 2*k+1);
            for al = k:L
                for bl = k:L-al
                    L_ma(al-k+1,bl-k+1) = al + bl;
                end
            end
            
            var   = 2.^(-beta*L_ma);
            var(1,1) = 1e-3;
            cost  = L_ma.*2.^L_ma;
            
            K_l = sum(sum( sqrt(var.*cost) ));
            M = max(2*ceil( eps^(-2)*K_l.*sqrt(var./cost) ),2);
            M = fliplr(triu(fliplr(M)));
            costl = 0;
            parfor j = 1:rel
                sumlsn  = 0;
                sumlret = 0;
                sumlreb = 0;
                % total degree index set
                for lx = k:L
                    for ly = k:L-lx
                        N = M(lx-k+1,ly-k+1);
                        %                 fprintf(1,'lx = %.1d, ly = %.1d \n',lx,ly)
                        [sums, cst] = mismc_l(lx,ly,N,data,rate,a,b,mu,k);
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
        sn = resn;
        re = rere;
        
        save('./workingdata/micomplexity_td.mat','sn','re','costs','L','k')
    case 2
        for i = 1:length(Eps)
            eps = Eps(i);
            alpha = 0.8;
            beta  = 1.6;

            L = ceil( log2(2/eps)/alpha );
            for al = k:L
                for bl = k:L
                    L_ma(al-k+1,bl-k+1) = al + bl;
                end
            end
            
            var   = 2.^(-beta*L_ma);
            var(1,1) = 1e-3;
            cost  = 2.^L_ma;
            
            K_l = sum(sum( sqrt(var.*cost) ));
            M = max(ceil(  eps^(-2)*K_l.*sqrt(var./cost) ),2);
            costl = 0;
            parfor j = 1:rel
                sumlsn  = 0;
                sumlret = 0;
                sumlreb = 0;
                % total degree index set
                for lx = k:L
                    for ly = k:L
                        N = M(lx-k+1,ly-k+1);
                        [sums, cst] = mismc_l(lx,ly,N,data,rate,a,b,mu,k);
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
        sn = resn;
        re = rere;
        
        save('./workingdata/micomplexity_tp.mat','sn','re','costs','L','k')        
end

end

