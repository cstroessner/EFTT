% This script generates the plot in fiugre 2

clear
load("data/GenzFuncs.mat")

figure

t = tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');

for ii = [1,3,5]
    nexttile
    semilogy(dlist1,(tabledata(ii,:,1)))
    hold on
    semilogy(dlist1,(tabledata(ii,:,2)))
    title(funcname{ii},'Interpreter','latex')
    xlabel('dimension $d$','Interpreter','latex')
    ylabel('$L^2$-error','Interpreter','latex')
    if ii == 1
        leg = legend('EFTT','DirectTT');
        set(leg,'Interpreter','latex','Location','northwest');
    end

    nexttile
    plot(dlist1,(tabledata(ii,:,3)))
    hold on
    plot(dlist1,(tabledata(ii,:,4)))
    title(funcname{ii},'Interpreter','latex')
    xlabel('dimension $d$','Interpreter','latex')
    ylabel('\# evals','Interpreter','latex')
end
set(gcf,'Position', [0 0 800 600])
print -depsc 'figure/Figure2'