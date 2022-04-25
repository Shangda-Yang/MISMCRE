Eps = 2.^(-8:1:-5);

alpha = 0.8;
temp = 0.5/log(2)/alpha;

logEps = -log(eps);

Levels = ceil(temp*(logEps + log(temp*logEps) + ...
    4*exp(log(2)*alpha*sqrt(2))*(temp^2+temp) ));

k = 5;
for lx = k:12
    for ly = k:12-lx
        fprintf(1,'Lx = %.1d and Ly = %.1d \n',lx,ly);
    end
end
