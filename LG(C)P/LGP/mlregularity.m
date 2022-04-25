close all; clear all;

randn('seed',20)
rand('seed',20)

s2=1.91;
b=(33/2/pi)^2;
% a=s2*gamma(3/2)/pi^(3/2)*sqrt(b);
% mu=log(126)-s2/2;
% b = 1;
a = 1;
mu = 0;

N    = 10;
M    = 20;
Lmax = 8;

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

rate = 2.6;

k = 5;
tic;
[we,se,cost] = ml_LGP_emp(@opre_mlLGP_l,N,M,Lmax,data,k,rate,a,b,mu);
toc;

save('./workingdata/ml_LGP_emperical.mat','Lmax','k',...
    'cost','we','se')

figure(1)
set(gcf,'Position',[300 300 1200 500]);
subplot(1,2,2)
plot(5:Lmax,log2(se(3,:)),'--or','Linewidth',2)
hold on
line = refline(-(rate-1),10);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(5:Lmax)
xlabel('$\ell_1,\ell_2$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 E[(\Delta F_{l,l}^{N} - \Delta F_{l,l})^2]$','Interpreter','latex','Fontsize',20);
legend('$l_1 = l_2$','$\mathcal{O}(2^{-1.6l_1})$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)

subplot(1,2,1)
plot(5:Lmax,log2(we(2,:)),'--or','Linewidth',2)
hold on
line = refline(-(rate-1)/2,0);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(5:Lmax)
xlabel('$\ell_1,\ell_2$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 |E[F_{l,l}-F_{l-1,l-1}]|$','Interpreter','latex','Fontsize',20);
legend('$l_1 = l_2$','$\mathcal{O}(2^{-0.8l_1})$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
saveas(gcf,'LGCMLSMC.fig')
saveas(gcf,'LGCMLSMC.jpg')