load('mlsmc_complexity.mat')

figure(1);
loglog(costs,mse_sn2,'linewidth',2)
hold on
loglog(costs,mse_re2,'linewidth',2)
loglog([1e5 1e7],[5e-5 5e-7],'k--')
grid on
xlabel('$Costs$','Interpreter','latex','Fontsize',15); 
ylabel('$log(MSE)$','Interpreter','latex','Fontsize',15);
legend('$SN$','$RE$','reference line with slop -1',...
    'Location','NorthEast','Interpreter','latex','Fontsize',15)
title('MLSMC','Interpreter','latex','Fontsize',15)

saveas(1,'MLSMCcomp.fig')
saveas(1,'MLSMCcomp.jpg')