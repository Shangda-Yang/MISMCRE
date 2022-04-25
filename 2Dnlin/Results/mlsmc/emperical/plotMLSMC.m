load('mlmc_emperical.mat')
set(gcf,'Position',[300 300 1200 500]);
% subplot(2,2,1)
% plot(k:L+k,log2(se1),'linewidth',2)
% hold on
% line = refline(-4, -2);
% line.LineStyle = '--';
% line.Color = 'k';
% grid on
% xticks(k:L+k)
% xlabel('$l$','Interpreter','latex','Fontsize',15); 
% ylabel('$log_2 E[(\Delta P_{l,l}^{N} - \Delta P_{l,l})^2]$','Interpreter','latex','Fontsize',15);
% legend('MLSMC','reference line with slop -4',...
%     'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Normalised Variance','Interpreter','latex','Fontsize',20)
% 
% subplot(2,2,2)
% plot(k:L+k,log2(we1),'linewidth',2)
% hold on
% line = refline(-2, -3);
% line.LineStyle = '--';
% line.Color = 'k';
% grid on
% xticks(k:L+k)
% xlabel('$l$','Interpreter','latex','Fontsize',15); 
% ylabel('$log_2 |E[P_{l,l}-P_{l-1,l-1}]|$','Interpreter','latex','Fontsize',15);
% legend('MLSMC','reference line with slop -2',...
%     'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Normalised Mean','Interpreter','latex','Fontsize',20)

subplot(1,2,2)
plot(k:L+k-1,log2(se2(1:end-1)),'--or','linewidth',2)
hold on
line = refline(-4, -5);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:L+k)
xlabel('$l$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 V_{l,l}$','Interpreter','latex','Fontsize',20);
legend(line,'$\mathcal{O}(-4l)$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Unnormalised Variance','Interpreter','latex','Fontsize',20)

subplot(1,2,1)
plot(k:L+k-1,log2(we2(1:end-1)),'--or','linewidth',2)
hold on
line = refline(-2, -3);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:L+k)
xlabel('$l$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 B_{l,l}$','Interpreter','latex','Fontsize',20);
legend(line,'$\mathcal{O}(-2l)$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Unnormalised Mean','Interpreter','latex','Fontsize',20)

saveas(gcf,'2DnlinMLSMC.fig')
saveas(gcf,'2DnlinMLSMC.jpg')