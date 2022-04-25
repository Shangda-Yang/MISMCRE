% ------------------------------------------------------ %
% function [u2] = analytic_u2(sigma,data,x_data)
% inputs:   sigma  = std deviation of the error of y - G
%           date   = observations
%           x_data = points of observations
% outputs:  u2     = analytic solution of u^2
% ------------------------------------------------------ %
function [u2] = analytic_u2(sigma,data,x_data)
x = x_data;
temp1 = sum(data.^2);
temp2 = sum(data.*(x.^2-x)/2);
temp3 = sum((x.^2-x).^2/4);
fun1 = @(z) exp(-(temp1 + 2*temp2.*z + temp3.*z.^2)/2/sigma/sigma);
cont = integral(fun1,-1, 1,'AbsTol',1e-16);
fun2 = @(z) z.^2.*exp(-(temp1 + 2*temp2.*z + temp3.*z.^2)/2/sigma/sigma)/cont;
u2 = integral(fun2,-1, 1,'AbsTol',1e-16);
end