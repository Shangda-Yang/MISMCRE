% ---------------------------------------------------------%
% function A = Multinomial_Resampling(w)
% Multinomial_Resampling is used for resampling
% input:  W = weights of samples
% output: A = sample indices with categorical distribution
% ---------------------------------------------------------%
function A = Multinomial_Resampling(w)
%     
%     N = length(W);
%     
%     s = W(1);
%     m = 1;
%     A = zeros(1,N);
%     
%     U = rand(1,N);
%     U = sort(U);
%     
%     for  n = 1:N
%         while s < U(n)
%             m = m + 1;
%             s = s + W(m);
%         end
%         A(n) = m;
%     end
    
M = length(w);
Q = cumsum(w);
Q(M)=1; % Just in case...
i=1;
while (i<=M)
    sampl = rand(1);
    j=1;
    while (Q(j)<sampl)
        j=j+1;
    end
    A(i)=j;
    i=i+1;
end
    
end