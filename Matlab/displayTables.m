% Display the resuts from generateData as tables


clear;clc;
%% Print results for all experiments except cookie
fprintf('\n\n-------------------------SUMMARY OF RESULTS-------------------------\n\n\n');

load('data/20BenchFuncsAdapt')
fprintf('-----------Table: EFTT (vs FTcross) for 20 benchmark functions with adaptive n-----------\n')
fprintf('|   function  |   error  |   eval  |   dof  | degree | TTrank | Tuckerrank |\n')
for ii = 1:n_funcs
    fprintf('| %11s | %8.2e | %7.i | %6.i |   %3.i  |   %2.i   |     %2.i     |\n',funcname{ii},...
        tabledata(ii,1),round(tabledata(ii,2)),round(tabledata(ii,3)),round(tabledata(ii,4)),...
        round(tabledata(ii,5)), round(tabledata(ii,6)))
end

clear
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

clear
load('data/20BenchFuncs')
fprintf('-----------Table: EFTT vs DirectTT for 20 benchmark functions with fixed n = 100-----------\n')
fprintf('|   function  |  method  |   error  |   eval  |   dof  | TTrank | Tuckerrank | reduce eval(%%) | reduce storage(%%) |\n')
for ii = 1:n_funcs
    fprintf('| %11s |   EFTT   | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |     %6.1f     |       %6.1f      |\n',funcname{ii},...
        tabledata(2*ii-1,1),round(tabledata(2*ii-1,2)),round(tabledata(2*ii-1,3)),round(tabledata(2*ii-1,4)),...
        round(tabledata(2*ii-1,5)),tablecompare(ii,1)*100,tablecompare(ii,2)*100)
    fprintf('|             | DirectTT | %8.2e | %7.i | %6.i |   %2.i   |            |                |                   |\n',...
        tabledata(2*ii,1),round(tabledata(2*ii,2)),round(tabledata(2*ii,3)),round(tabledata(2*ii,4)))
end

clear
load('data/GenzFuncs')
fprintf('-----------Table: EFTT vs DirectTT for Genz functions with fixed n = 100-----------\n')
for ii = 1:n_funcs
    fprintf('|    function   |  d  |  method  |   error  |   eval  |   dof  | TTrank | Tuckerrank | reduce eval(%%) | reduce storage(%%) |\n')
    for jj = 1:length(dlist1)
        d = dlist1(jj);
        if ii == 2
            d = dlist2(jj);
        end

        if jj == 4
            fprintf('| %13s | %3.i |   EFTT   | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |     %6.1f     |       %6.1f      |\n',...
                funcname{ii},d,tabledata(ii,jj,1),round(tabledata(ii,jj,3)),round(tabledata(ii,jj,5)),round(tabledata(ii,jj,7)),...
                round(tabledata(ii,jj,9)),tablecompare(ii,jj,1)*100,tablecompare(ii,jj,2)*100)
        else
            fprintf('|               | %3.i |   EFTT   | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |     %6.1f     |       %6.1f      |\n',...
                d,tabledata(ii,jj,1),round(tabledata(ii,jj,3)),round(tabledata(ii,jj,5)),round(tabledata(ii,jj,7)),...
                round(tabledata(ii,jj,9)),tablecompare(ii,jj,1)*100,tablecompare(ii,jj,2)*100)
        end

        fprintf('|               |     |   FTT    | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |                |                   |\n',...
            tabledata(ii,jj,2),round(tabledata(ii,jj,4)),round(tabledata(ii,jj,6)),round(tabledata(ii,jj,8)),...
            round(tabledata(ii,jj,10)))
    end
    fprintf('------------------------------------------------------------------------\n')
end