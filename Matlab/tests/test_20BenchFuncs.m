function [tabledata,tablecompare] = test_20BenchFuncs(n_repetitions,savedata)
% Test 20 benchmark functions for DirectTT and EFTT with fixed n = 100
% The results in the paper are obtained by running the experiment 100 times
% and taking an average. You may simply set n_repetitions to 1 to get the
% table quickly.

%savedata = false;
%n_repetitions = 100;
tol = 1e-10;
maxTTranks = 50;
maxTuckerranks = 50;

n_funcs = 20;
dlist = [7*ones(1,12),8,6,8,10,5,6,8,3];
funcname = {'ackley','alpine','dixon','exponential','grienwank',...
    'michalewicz','piston','qing','rastrigin','rosenbrock',...
    'schaffer','schwefel','borehole','otlcircuit','robot',...
    'wingweight','fried','grlee09','detpep108d','detpep10exp'};


load("testPoints.mat",'x_all','dlist_all')
saveExceptThem = '^(?!(x_all|dlist_all)$).';

for kk = 1:n_repetitions
    rng(kk)
    fprintf('-----------repetition = %i-----------\n',kk)
    fprintf('|   function  |  method  |   error  |   eval  |   dof  | TTrank | Tuckerrank | reduce eval(%%) | reduce storage(%%) |\n')
    for ii = 1:n_funcs
        d = dlist(ii);

        
        x = x_all{dlist_all == d};

        selectedFunction = ii;
        y = applyToRows(@(x)testfunction(x,'selectedfunction',selectedFunction),x);

        n = 100*ones(d,1);
        eftt = EFTT(@(x)testfunction(x,'selectedfunction',selectedFunction),d,'n',n,'tol',tol,'ifadapt',zeros(d,1),...
            'TuckerRanks',maxTuckerranks,'ttranks',maxTTranks,'verb',0);
        EFTTerror(kk,ii) = norm(fevaluate(eftt,x) - y)/norm(y);
        EFTTeval(kk,ii) = eftt.numAllEvals;
        EFTTdof(kk,ii) = dof(eftt);
        EFTTttranks{kk,ii} = ttrank(eftt);
        EFTTtuckerranks{kk,ii} = tuckerrank(eftt);


        ftt = FTT(@(x)testfunction(x,'selectedfunction',selectedFunction),d,'n',n,'tol',tol,'ttranks',maxTTranks,'verb',0);
        FTTerror(kk,ii) = norm(fevaluate(ftt,x) - y)/norm(y);
        FTTeval(kk,ii) = ftt.numAllEvals;
        FTTdof(kk,ii) = dof(ftt);
        FTTttranks{kk,ii} = ttrank(ftt);

        CompareEval(kk,ii) = (1-EFTTeval(kk,ii)/FTTeval(kk,ii))*100;
        CompareDof(kk,ii) = (1-EFTTdof(kk,ii)/FTTdof(kk,ii))*100;

        fprintf('| %11s |   EFTT   | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |     %6.1f     |       %6.1f      |\n',funcname{ii},...
            EFTTerror(kk,ii),EFTTeval(kk,ii),EFTTdof(kk,ii),max(EFTTttranks{kk,ii}),...
            max(EFTTtuckerranks{kk,ii}),CompareEval(kk,ii),CompareDof(kk,ii))
        fprintf('|             | DirectTT | %8.2e | %7.i | %6.i |   %2.i   |            |                |                   |\n',...
            FTTerror(kk,ii),FTTeval(kk,ii),FTTdof(kk,ii),max(FTTttranks{kk,ii}))
    end
    if savedata
        save('data/20BenchFuncs','-regexp', saveExceptThem)
    end
end


%error, eval, dof, ttranks, tuckerranks
tabledata = zeros(2*n_funcs,5);
for ii = 1:n_funcs
    tabledata(2*ii-1,1) = geomean(EFTTerror(:,ii));
    tabledata(2*ii-1,2) = mean(EFTTeval(:,ii));
    tabledata(2*ii-1,3) = mean(EFTTdof(:,ii));
    tempttrank = cell2mat(EFTTttranks(:,ii)');
    tabledata(2*ii-1,4) = max(mean(tempttrank(2:end-1,:),1));
    tabledata(2*ii-1,5) = max(mean(cell2mat(EFTTtuckerranks(:,ii)'),1));

    tabledata(2*ii,1) = geomean(FTTerror(:,ii));
    tabledata(2*ii,2) = mean(FTTeval(:,ii));
    tabledata(2*ii,3) = mean(FTTdof(:,ii));
    tempttrank = cell2mat(FTTttranks(:,ii)');
    tabledata(2*ii,4) = max(mean(tempttrank,1));
    tabledata(2*ii,5) = 0;
end
% reduce eval (%), reduce storage (%), increse accuracy (digit)
tablecompare = zeros(n_funcs,3);
tablecompare(:,1) = 1 - tabledata(1:2:2*n_funcs,2)./tabledata(2:2:2*n_funcs,2);
tablecompare(:,2) = 1 - tabledata(1:2:2*n_funcs,3)./tabledata(2:2:2*n_funcs,3);
tablecompare(:,3) = log10(tabledata(2:2:2*n_funcs,4)./tabledata(1:2:2*n_funcs,4));

if savedata
    save('data/20BenchFuncs','-regexp', saveExceptThem)
end
%load('data/20BenchFuncs')
fprintf('-----------Table: EFTT vs DirectTT for 20 benchmark functions with fixed n = 100-----------\n')
fprintf('|   function  |  method  |   error  |   eval  |   dof  | TTrank | Tuckerrank | reduce eval(%%) | reduce storage(%%) |\n')
for ii = 1:n_funcs
    fprintf('| %11s |   EFTT   | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |     %6.1f     |       %6.1f      |\n',funcname{ii},...
        tabledata(2*ii-1,1),round(tabledata(2*ii-1,2)),round(tabledata(2*ii-1,3)),round(tabledata(2*ii-1,4)),...
        round(tabledata(2*ii-1,5)),tablecompare(ii,1)*100,tablecompare(ii,2)*100)
    fprintf('|             | DirectTT | %8.2e | %7.i | %6.i |   %2.i   |            |                |                   |\n',...
        tabledata(2*ii,1),round(tabledata(2*ii,2)),round(tabledata(2*ii,3)),round(tabledata(2*ii,4)))
end
