load('mismc_complexity_tp.mat')
mse_sn_tp = mse_sn2;
mse_re_tp = mse_re2;
costs_tp  = costs;

load('mismc_complexity_td.mat')
mse_sn_td = mse_sn2;
mse_re_td = mse_re2;
costs_td  = costs;

figure(1);
set(gcf,'Position',[300 300 1200 1000]);
subplot(2,2,1)
loglog(costs_tp,mse_sn_tp,'linewidth',2)
hold on
loglog(costs_tp,mse_re_tp,'linewidth',2)
loglog([1e4 1e8],[5e-4 5e-8],'k--')
grid on
xlabel('$Costs$','Interpreter','latex','Fontsize',15); 
ylabel('$MSE$','Interpreter','latex','Fontsize',15);
legend('$SN_{tp}$','$RE_{tp}$','reference line with slop -1',...
    'Location','SouthWest','Interpreter','latex','Fontsize',15)


subplot(2,2,2)
loglog(costs_td,mse_sn_td,'linewidth',2)
hold on
loglog(costs_td,mse_re_td,'linewidth',2)
loglog([1e4 1e8],[5e-4 5e-8],'k--')
grid on
xlabel('$Costs$','Interpreter','latex','Fontsize',15); 
ylabel('$MSE$','Interpreter','latex','Fontsize',15);
legend('$SN_{td}$','$RE_{td}$','reference line with slop -1',...
    'Location','SouthWest','Interpreter','latex','Fontsize',15)

subplot(2,2,3)
loglog(costs_tp,mse_sn_tp,'linewidth',2)
hold on
loglog(costs_td,mse_sn_td,'linewidth',2)
loglog([1e4 1e8],[5e-4 5e-8],'k--')
grid on
xlabel('$Costs$','Interpreter','latex','Fontsize',15); 
ylabel('$MSE$','Interpreter','latex','Fontsize',15);
legend('$SN_{tp}$','$SN_{td}$','reference line with slop -1',...
    'Location','SouthWest','Interpreter','latex','Fontsize',15)

subplot(2,2,4)
loglog(costs_tp,mse_re_tp,'linewidth',2)
hold on
loglog(costs_td,mse_re_td,'linewidth',2)
loglog([1e4 1e8],[5e-4 5e-8],'k--')
grid on
xlabel('$Costs$','Interpreter','latex','Fontsize',15); 
ylabel('$MSE$','Interpreter','latex','Fontsize',15);
legend('$RE_{tp}$','$RE_{td}$','reference line with slop -1',...
    'Location','SouthWest','Interpreter','latex','Fontsize',15)

sgtitle('MISMC','Interpreter','latex','Fontsize',20)

saveas(1,'MISMCcomp.fig')
saveas(1,'MISMCcomp.jpg')
