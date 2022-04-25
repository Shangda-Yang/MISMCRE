close all; clear all;

% nslots = str2double(getenv('NSLOTS'));
% parpool(nslots);
% pp = gcp;
% poolsize = pp.NumWorkers;

randn('seed',20)
rand('seed',20)

% s2:   \sigma^2
% mu:   \theta_1
% a:    \theta_2
% b:    \theta_3
% rate: \beta^{\prime}

s2=1.91;
b=(33/2/pi)^2;
% a=s2*gamma(3/2)/pi^(3/2)*sqrt(b);
% mu=log(126)-s2/2;
% b = 1;
a = 1;
mu = 0;

N    = 100; % number of samples
M    = 20;   % number of realisations
Lmax = 8;    % Lmax maximum level for emperical tests

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

rate = 2.6;

k = 5;
tic;
[we,se,cost] = mi_LGP_emp(@opre_miLGP_l,N,M,Lmax,data,k,rate,a,b,mu);
toc;

% save('./workingdata/mi_LGP_emperical.mat','Lmax','k',...
%     'cost','we','se')

test1 = se{2};
test2 = we{2};

subplot(3,2,1)
contourf(k+1:Lmax,k+1:Lmax,test1(2:end,2:end))
xticks(k+1:Lmax)
yticks(k+1:Lmax)
xlabel('$lx$','Interpreter','latex','Fontsize',20); 
ylabel('$ly$','Interpreter','latex','Fontsize',20); 
colorbar

subplot(3,2,2)
contourf(k+1:Lmax,k+1:Lmax,test2(2:end,2:end))
xticks(k+1:Lmax)
yticks(k+1:Lmax)
xlabel('$lx$','Interpreter','latex','Fontsize',20); 
ylabel('$ly$','Interpreter','latex','Fontsize',20); 
colorbar

subplot(3,2,3)
plot(k:Lmax,log2(diag(test1)),'Linewidth',2)
hold on
line = refline(-2*(rate-1),0);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lmax)
xlabel('$l$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 MSE$','Interpreter','latex','Fontsize',20);
legend('MSE','-3.2L',...
    'Location','NorthEast','Interpreter','latex','Fontsize',20)
title('Regularity','Interpreter','latex','Fontsize',20)

subplot(3,2,4)
plot(k:Lmax,log2(diag(test2)),'Linewidth',2)
hold on
line = refline(-(rate-1),0);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lmax)
xlabel('$l$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 MSE$','Interpreter','latex','Fontsize',20);
legend('MSE','-1.6L',...
    'Location','NorthEast','Interpreter','latex','Fontsize',20)
title('Regularity','Interpreter','latex','Fontsize',20)

subplot(3,2,5)
plot(k:Lmax,log2(test1(end,:)),'Linewidth',2)
hold on
line = refline(-(rate-1),0);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lmax)
xlabel('$l$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 MSE$','Interpreter','latex','Fontsize',20);
legend('MSE','-1.6L',...
    'Location','NorthEast','Interpreter','latex','Fontsize',20)
title('Regularity','Interpreter','latex','Fontsize',20)

subplot(3,2,6)
plot(k:Lmax,log2(test2(end,:)),'Linewidth',2)
hold on
line = refline(-(rate-1)/2,0);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lmax)
xlabel('$l$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 MSE$','Interpreter','latex','Fontsize',20);
legend('MSE','-0.8L',...
    'Location','NorthEast','Interpreter','latex','Fontsize',20)
title('Regularity','Interpreter','latex','Fontsize',20)
