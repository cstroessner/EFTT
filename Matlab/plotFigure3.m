
clear
close all
load('data/SineIntegrationAdapt')
fprintf('-----------Table: EFTT (vs FTcross) for sine function integration with adaptive n-----------\n')
fprintf('|  d  |   tol   |   EFTT   |   error  |  eval |  dof  | degree |\n')
for ii = 1:length(dlist)
    d = dlist(ii);
    for jj = 1:length(tollist)
        tol = tollist(jj);
        if jj == 2
            fprintf('| %3.i | %7.1e |   EFTT   | %8.2e | %5.i | %5.i |   %2.i   |\n',d,tol,...
                tabledata(ii,jj,1),round(tabledata(ii,jj,2)),round(tabledata(ii,jj,3)),round(tabledata(ii,jj,4)));
        else
            fprintf('|     | %7.1e |   EFTT   | %8.2e | %5.i | %5.i |   %2.i   |\n',tol,...
                tabledata(ii,jj,1),round(tabledata(ii,jj,2)),round(tabledata(ii,jj,3)),round(tabledata(ii,jj,4)));
        end
    end
    fprintf('------------------------------------------------------------------\n')
end

% values for c3py generated in python and copied directly in the plots below


%% plot 1 error vs dim

close all

figure(8)
set(gca,'fontsize',10)
set(figure(8), 'Position', [0 0 370 300])
semilogy(dlist,[2.16780106e-14 1.53708369e-13 1.49273653e-12 1.08012368e-12  1.13350993e-12 1.16476755e-12 1.07746807e-12 2.28356113e-12 2.44428866e-12 2.58998307e-12],'r') %low tol
hold on
semilogy(dlist,tabledata(:,1,1),'b')
xlabel('dimension $d$','Interpreter','latex')
ylabel('error','Interpreter','latex')
leg = legend('c3py','EFTT');
set(leg,'Interpreter','latex','location','northwest');
print -depsc 'figure/SinIntPlot1'


% plot evals vs dim
figure(9)
set(gca,'fontsize',10)
set(figure(9), 'Position', [0 0 370 300])
plot(dlist,[ 1512.  9618. 17934. 26292. 34692. 42945. 51464. 59892. 68446. 76671.],'r') %low tol
hold on
plot(dlist,tabledata(:,1,2),'b')
xlabel('dimension $d$','Interpreter','latex')
ylabel('\# evals','Interpreter','latex')
leg = legend('c3py','EFTT');
set(leg,'Interpreter','latex','location','northwest');
print -depsc 'figure/SinIntPlot2'

% plot dofs vs dim
figure(10)
set(gca,'fontsize',10)
set(figure(10), 'Position', [0 0 370 300])
plot(dlist,[  417.  2625.  4879.  7133.  9387. 11641. 13895. 16149. 18403. 20657.],'r') %low tol
hold on
plot(dlist,tabledata(:,1,3),'b')
xlabel('dimension $d$','Interpreter','latex')
ylabel('\# dofs','Interpreter','latex')
leg = legend('c3py','EFTT');
set(leg,'Interpreter','latex','location','northwest');
print -depsc 'figure/SinIntPlot2'
