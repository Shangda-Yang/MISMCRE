% ------------------------------------------------------ %
% function [u_vector] = FEM_u(l,uni)
% finite-element u_vector of Au = f
% inputs:   l         = level of refinement
%           uni       = random input   
% outputs:  u_vector  = FEM grid solutions
% ------------------------------------------------------ %
function [u_vector] = FEM_u(l,uni)    
    
    h = 2^-l;
    
    K = 1/h;
    
    invh = h^(-1);
    
    A_i_j   = 2*ones(K-1,1); % diagonal
    A_i_jm1 = -ones(K-1,1); % lower
    A_i_jp1 = -ones(K-1,1); % upper
    
    A = spdiags([A_i_jm1 A_i_j A_i_jp1],[-1 0 1],K-1,K-1);
    
    A = A.*invh;
    
    f = uni*h*ones(K-1,1);
    
    u_vector = [0; A\f; 0];
end
       
    
    