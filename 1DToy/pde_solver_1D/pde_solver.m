% ------------------------------------------------------ %
% function [v]  = pde_solver(x,l,uni)
% inputs:   x   = interpolation points
%           l   = level of refinement
%           uni = random input
% outputs:  v   = solution at corresponding points
% ------------------------------------------------------ %

function [v] = pde_solver(x,l,uni)
if l == 0
   
    v = 0;

else
    h = 2^-l;
    u_vector = FEM_u(l,uni);
    x_vector = (0:h:1);
    v = zeros(length(x),1);
    for i = 1:length(x)
        if x(i) == 1 || x(i) == 0
            v(i) = 0;
        else
            ind = find(x(i) <= x_vector, 1);
            x_u = x_vector(ind);
            x_d = x_vector(ind - 1);
        
            v(i) = ( x_u - x(i) )/h*u_vector(ind-1) + ...
               ( x(i) - x_d )/h*u_vector(ind);
        end
    end
end
end
