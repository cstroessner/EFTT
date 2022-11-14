% This script generates the data for the tables in the appendix and the raw
% data for Figures 2 and 3

clear;
rng(1)

%% User specified variables
savedata = true;
n_realizations = 30;
n_repetitions = 100;

run_experiments = {'20BenchFuncsAdapt','20BenchFuncs','SineIntegrationAdapt','GrenzFuncs'};

%% run experiments

diary diary_test_experiments

if ismember('20BenchFuncsAdapt',run_experiments)
    fprintf('----------Running test_20BenchFuncsAdapt----------\n');
    tabledata{1,1} = test_20BenchFuncsAdapt(1,savedata);
end


if ismember('SineIntegrationAdapt',run_experiments)
    fprintf('----------Running test_SineIntegrationAdapt----------\n');
    tabledata{3,1} = test_SineIntegrationAdapt(n_repetitions,savedata);
end

if ismember('20BenchFuncs',run_experiments)
    fprintf('----------Running test_20BenchFuncs----------\n');
    [tabledata{4,1},tablecompare{4,1}] = test_20BenchFuncs(n_repetitions,savedata);
end

if ismember('GenzFuncs',run_experiments)
    fprintf('----------Running test_GenzFuncs----------\n');
    [tabledata{5,1},tablecompare{5,1}] = test_GenzFuncs(n_realizations,savedata);
end


%% Print results
fprintf('\n\n-------------------------SUMMARY OF RESULTS-------------------------\n\n\n');

if ismember('20BenchFuncsAdapt',run_experiments)
    n_funcs = 20;
    funcname = {'ackley','alpine','dixon','exponential','grienwank',...
        'michalewicz','piston','qing','rastrigin','rosenbrock',...
        'schaffer','schwefel','borehole','otlcircuit','robot',...
        'wingweight','fried','grlee09','detpep108d','detpep10exp'};

    fprintf('-----------Table: EFTT (vs FTcross) for 20 benchmark functions with adaptive n-----------\n')
    fprintf('|   function  |   error  |   eval  |   dof  | degree | TTrank | Tuckerrank |\n')
    for ii = 1:n_funcs
        fprintf('| %11s | %8.2e | %7.i | %6.i |   %3.i  |   %2.i   |     %2.i     |\n',funcname{ii},...
            tabledata{1,1}(ii,1),round(tabledata{1,1}(ii,2)),round(tabledata{1,1}(ii,3)),round(tabledata{1,1}(ii,4)),...
            round(tabledata{1,1}(ii,5)), round(tabledata{1,1}(ii,6)))
    end
end

if ismember('SineIntegrationAdapt',run_experiments)
    dlist = linspace(10,685,10);
    tollist = logspace(-13, -1, 10);
    tollist = tollist([1,5,10]);

    fprintf('-----------Table: EFTT (vs FTcross) for sine function integration with adaptive n-----------\n')
    fprintf('|  d  |   tol   |   EFTT   |   error  |  eval |  dof  | degree |\n')
    for ii = 1:length(dlist)
        d = dlist(ii);
        for jj = 1:length(tollist)
            tol = tollist(jj);
            if jj == 2
                fprintf('| %3.i | %7.1e |   EFTT   | %8.2e | %5.i | %5.i |   %2.i   |\n',d,tol,...
                    tabledata{3,1}(ii,jj,1),round(tabledata{3,1}(ii,jj,2)),round(tabledata{3,1}(ii,jj,3)),round(tabledata{3,1}(ii,jj,4)));
            else
                fprintf('|     | %7.1e |   EFTT   | %8.2e | %5.i | %5.i |   %2.i   |\n',tol,...
                    tabledata{3,1}(ii,jj,1),round(tabledata{3,1}(ii,jj,2)),round(tabledata{3,1}(ii,jj,3)),round(tabledata{3,1}(ii,jj,4)));
            end
        end
        fprintf('------------------------------------------------------------------\n')
    end
end

if ismember('20BenchFuncs',run_experiments)
    n_funcs = 20;
    funcname = {'ackley','alpine','dixon','exponential','grienwank',...
        'michalewicz','piston','qing','rastrigin','rosenbrock',...
        'schaffer','schwefel','borehole','otlcircuit','robot',...
        'wingweight','fried','grlee09','detpep108d','detpep10exp'};

    fprintf('-----------Table: EFTT vs DirectTT for 20 benchmark functions with fixed n = 100-----------\n')
    fprintf('|   function  |  method  |   error  |   eval  |   dof  | TTrank | Tuckerrank | reduce eval(%%) | reduce storage(%%) |\n')
    for ii = 1:n_funcs
        fprintf('| %11s |   EFTT   | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |     %6.1f     |       %6.1f      |\n',funcname{ii},...
            tabledata{4,1}(2*ii-1,1),round(tabledata{4,1}(2*ii-1,2)),round(tabledata{4,1}(2*ii-1,3)),round(tabledata{4,1}(2*ii-1,4)),...
            round(tabledata{4,1}(2*ii-1,5)),tablecompare{4,1}(ii,1)*100,tablecompare{4,1}(ii,2)*100)
        fprintf('|             | DirectTT | %8.2e | %7.i | %6.i |   %2.i   |            |                |                   |\n',...
            tabledata{4,1}(2*ii,1),round(tabledata{4,1}(2*ii,2)),round(tabledata{4,1}(2*ii,3)),round(tabledata{4,1}(2*ii,4)))
    end
end

if ismember('GenzFuncs',run_experiments)
    dlist1 = [20,50,100,200,300,400,500];
    dlist2 = [5,10,15,20,25,30,50];
    funcname = {'Oscillatory','ProductPeak','CornerPeak','Gaussian','Continuous','Discontinuous'};
    n_funcs = 6;

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
                    funcname{ii},d,tabledata{5,1}(ii,jj,1),round(tabledata{5,1}(ii,jj,3)),round(tabledata{5,1}(ii,jj,5)),round(tabledata{5,1}(ii,jj,7)),...
                    round(tabledata{5,1}(ii,jj,9)),tablecompare{5,1}(ii,jj,1)*100,tablecompare{5,1}(ii,jj,2)*100)
            else
                fprintf('|               | %3.i |   EFTT   | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |     %6.1f     |       %6.1f      |\n',...
                    d,tabledata{5,1}(ii,jj,1),round(tabledata{5,1}(ii,jj,3)),round(tabledata{5,1}(ii,jj,5)),round(tabledata{5,1}(ii,jj,7)),...
                    round(tabledata{5,1}(ii,jj,9)),tablecompare{5,1}(ii,jj,1)*100,tablecompare{5,1}(ii,jj,2)*100)
            end

            fprintf('|               |     |   FTT    | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |                |                   |\n',...
                tabledata{5,1}(ii,jj,2),round(tabledata{5,1}(ii,jj,4)),round(tabledata{5,1}(ii,jj,6)),round(tabledata{5,1}(ii,jj,8)),...
                round(tabledata{5,1}(ii,jj,10)))
        end
        fprintf('------------------------------------------------------------------------\n')
    end
end
diary off