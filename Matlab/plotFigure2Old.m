% This script generates the plot in fiugre 2

clear
load("data/GenzFuncs.mat")

%% error plot for function
close all

figure(8)
set(gca,'fontsize',10)
set(figure(8), 'Position', [0 0 370 300])
for i = [1,3,5]
semilogy(dlist1,(tabledata(i,:,1)./tabledata(i,:,2)))
hold on
end
xlabel('dimension $d$','Interpreter','latex')
ylabel('$L^2$-error EFTT / $L^2$-error DirectTT','Interpreter','latex')
leg = legend('oscillatory','corner peak','continuous');
set(leg,'Interpreter','latex','location','northwest');
print -depsc 'figure/RankImpactError'

%% eval plot for function
figure(7)
set(gca,'fontsize',10)
set(figure(7), 'Position', [0 0 370 300])
for i = [1,3,5]
plot(dlist1,tabledata(i,:,3)./tabledata(i,:,4))
hold on
end
xlabel('dimension $d$','Interpreter','latex')
ylabel('evals EFTT / evals DirectTT','Interpreter','latex')
leg = legend('oscillatory','corner peak','continuous');
ylim([0.1,1])
set(leg,'Interpreter','latex','Location','southeast');
print -depsc 'figure/RankImpactEvals'