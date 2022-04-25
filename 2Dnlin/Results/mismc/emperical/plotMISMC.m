load('mismc_emperical.mat')
%% lx = ly
figure(1);
set(gcf,'Position',[600 600 1200 500]);
% subplot(2,2,1)
% plot(k:Lx+k,log2(diag(se1)),'linewidth',2)
% hold on
% line = refline(-8, -2);
% line.LineStyle = '--';
% line.Color = 'k';
% grid on
% xticks(k:Lx+k)
% xlabel('$l$','Interpreter','latex','Fontsize',15); 
% ylabel('$log_2 E[(\Delta P_{l,l}^{N} - \Delta P_{l,l})^2]$','Interpreter','latex','Fontsize',15);
% legend('MISMC','reference line with slop -8',...
%     'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Normalised Variance','Interpreter','latex','Fontsize',20)

% subplot(2,2,2)
% plot(k:Lx+k,log2(diag(we1)),'linewidth',2)
% hold on
% line = refline(-4, -3);
% line.LineStyle = '--';
% line.Color = 'k';
% grid on
% xticks(k:Lx+k)
% xlabel('$l$','Interpreter','latex','Fontsize',15); 
% ylabel('$log_2 |E[\Delta P]|$','Interpreter','latex','Fontsize',15);
% legend('MISMC','reference line with slop -4',...
%     'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Normalised Mean','Interpreter','latex','Fontsize',20)

subplot(1,2,2)
plot(k:Lx+k-1,log2(diag(se2(1:end-1,1:end-1))),'--or','linewidth',2)
hold on
line = refline(-8, -3);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lx+k-1)
xlim([k Lx+k-1])
xlabel('$\alpha_1,\alpha_2$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 V_{\alpha_1,\alpha_2}$','Interpreter','latex','Fontsize',20);
legend('$\alpha_1=\alpha_2$','$-8\alpha_1$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Unnormalised Variance','Interpreter','latex','Fontsize',20)

subplot(1,2,1)
plot(k:Lx+k-1,log2(diag(we2(1:end-1,1:end-1))),'--or','linewidth',2)
hold on
line = refline(-4, -3);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lx+k-1)
xlim([k Lx+k-1])
xlabel('$\alpha_1,\alpha_2$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 B_{\alpha_1,\alpha_2}$','Interpreter','latex','Fontsize',20);
legend('$\alpha_1=\alpha_2$','$-4\alpha_1$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Unnormalised Mean','Interpreter','latex','Fontsize',20)

% sgtitle('lx = ly','Interpreter','latex','Fontsize',20)
saveas(1,'2DnlinMISMC12.fig')
saveas(1,'2DnlinMISMC12.jpg')
%% fixed lx
figure(2);
set(gcf,'Position',[600 600 1200 500]);
% subplot(2,2,1)
% plot(k:Ly+k,log2(se1(:,end)),'linewidth',2)
% hold on
% line = refline(-4, -20);
% line.LineStyle = '--';
% line.Color = 'k';
% grid on
% xticks(k:Lx+k)
% xlabel('$l$','Interpreter','latex','Fontsize',15); 
% ylabel('$log_2 E[(\Delta P_{l,l}^{N} - \Delta P_{l,l})^2]$','Interpreter','latex','Fontsize',15);
% legend('MISMC','reference line with slop -4',...
%     'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Normalised Variance','Interpreter','latex','Fontsize',20)

% subplot(2,2,2)
% plot(k:Ly+k,log2(we1(:,end)),'linewidth',2)
% hold on
% line = refline(-2, -15);
% line.LineStyle = '--';
% line.Color = 'k';
% grid on
% xticks(k:Lx+k)
% xlabel('$l$','Interpreter','latex','Fontsize',15); 
% ylabel('$log_2 |E[\Delta P]|$','Interpreter','latex','Fontsize',15);
% legend('MISMC','reference line with slop -2',...
%     'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Normalised Mean','Interpreter','latex','Fontsize',20)

subplot(1,2,2)
plot(k:Ly+k-1,log2(se2(1:end-1,end-1)),'--or','linewidth',2)
hold on
line = refline(-4, -30);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lx+k-1)
xlim([k Lx+k-1])
xlabel('$\alpha_2$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 V_{\alpha_1,\alpha_2}$','Interpreter','latex','Fontsize',20);
legend('$\alpha_1 = 6$','$-4\alpha_2$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Unnormalised Variance','Interpreter','latex','Fontsize',20)

subplot(1,2,1)
plot(k:Ly+k-1,log2(we2(1:end-1,end-1)),'--or','linewidth',2)
hold on
line = refline(-2, -20);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lx+k-1)
xlim([k Lx+k-1])
xlabel('$\alpha_2$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 B_{\alpha_1,\alpha_2}$','Interpreter','latex','Fontsize',20);
legend('$\alpha_1 = 6$','$-2\alpha_2$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Unnormalised Mean','Interpreter','latex','Fontsize',20)

% sgtitle('fixed lx','Interpreter','latex','Fontsize',20)
saveas(2,'2DnlinMISMC2.fig')
saveas(2,'2DnlinMISMC2.jpg')
%% fixed ly
figure(3);
set(gcf,'Position',[600 600 1200 500]);
% subplot(2,2,1)
% plot(k:Lx+k,log2(se1(end,:)),'linewidth',2)
% hold on
% line = refline(-4, -20);
% line.LineStyle = '--';
% line.Color = 'k';
% grid on
% xticks(k:Lx+k)
% xlabel('$l$','Interpreter','latex','Fontsize',15); 
% ylabel('$log_2 E[(\Delta P_{l,l}^{N} - \Delta P_{l,l})^2]$','Interpreter','latex','Fontsize',15);
% legend('MISMC','reference line with slop -4',...
%     'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Normalised Variance','Interpreter','latex','Fontsize',20)

% subplot(2,2,2)
% plot(k:Lx+k,log2(we1(end,:)),'linewidth',2)
% hold on
% line = refline(-2, -15);
% line.LineStyle = '--';
% line.Color = 'k';
% grid on
% xticks(k:Lx+k)
% xlabel('$l$','Interpreter','latex','Fontsize',15); 
% ylabel('$log_2 |E[\Delta P]|$','Interpreter','latex','Fontsize',15);
% legend('MISMC','reference line with slop -2',...
%     'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Normalised Mean','Interpreter','latex','Fontsize',20)

subplot(1,2,2)
plot(k:Lx+k-1,log2(se2(end-1,1:end-1)),'--or','linewidth',2)
hold on
line = refline(-4, -30);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lx+k-1)
xlim([k Lx+k-1])
xlabel('$\alpha_1$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 V_{\alpha_1,\alpha_2}$','Interpreter','latex','Fontsize',20);
legend('$\alpha_2 = 6$','$-4\alpha_1$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Unnormalised Variance','Interpreter','latex','Fontsize',20)

subplot(1,2,1)
plot(k:Lx+k-1,log2(we2(end-1,1:end-1)),'--or','linewidth',2)
hold on
line = refline(-2, -20);
line.LineStyle = '--';
line.Color = 'k';
grid on
xticks(k:Lx+k-1)
xlim([k Lx+k-1])
xlabel('$\alpha_1$','Interpreter','latex','Fontsize',20); 
ylabel('$log_2 B_{\alpha_1,\alpha_2}$','Interpreter','latex','Fontsize',20);
legend('$\alpha_2 = 6$','$-2\alpha_1$',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
% title('Unnormalised Mean','Interpreter','latex','Fontsize',20)
% 
% sgtitle('fixed ly','Interpreter','latex','Fontsize',20)
saveas(3,'2DnlinMISMC1.fig')
saveas(3,'2DnlinMISMC1.jpg')
%% contour plot
figure(4)
set(gcf,'Position',[600 600 1200 500]);
% subplot(2,2,1)
% contourf(k:Lx+k,k:Ly+k,log2(se1))
% xticks(k:Lx+k)
% yticks(k:Ly+k)
% xlabel('$lx$','Interpreter','latex','Fontsize',15); 
% ylabel('$ly$','Interpreter','latex','Fontsize',15); 
% title('Normalised Variance','Interpreter','latex','Fontsize',20)
% colorbar
% 
% subplot(2,2,2)
% contourf(k:Lx+k,k:Ly+k,log2(we1))
% xticks(k:Lx+k)
% yticks(k:Ly+k)
% xlabel('$l_x$','Interpreter','latex','Fontsize',15); 
% ylabel('$l_y$','Interpreter','latex','Fontsize',15); 
% title('Normalised Mean','Interpreter','latex','Fontsize',20)
% colorbar

subplot(1,2,2)
contourf(k:Lx+k-1,k:Ly+k-1,log2(se2(1:end-1,1:end-1)))
xticks(k:Lx+k-1)
yticks(k:Ly+k-1)
xlabel('$\alpha_1$','Interpreter','latex','Fontsize',20); 
ylabel('$\alpha_2$','Interpreter','latex','Fontsize',20); 
% title('Unnormalised Variance','Interpreter','latex','Fontsize',20)
colorbar

subplot(1,2,1)
contourf(k:Lx+k-1,k:Ly+k-1,log2(we2(1:end-1,1:end-1)))
xticks(k:Lx+k-1)
yticks(k:Ly+k-1)
xlabel('$\alpha_1$','Interpreter','latex','Fontsize',20); 
ylabel('$\alpha_2$','Interpreter','latex','Fontsize',20); 
% title('Unnormalised Mean','Interpreter','latex','Fontsize',20)
colorbar

saveas(4,'2DnlinContour.fig')
saveas(4,'2DnlinContour.jpg')