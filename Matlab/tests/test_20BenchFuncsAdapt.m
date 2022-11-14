function tabledata = test_20BenchFuncsAdapt(n_repetitions,savedata)
% Test 20 benchmark functions for EFTT with n choosen adaptively
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
    fprintf('|   function  |   error  |   eval  |   dof  | degree | TTrank | Tuckerrank |\n')
    for ii = 1:n_funcs
        d = dlist(ii);

        
        x = x_all{dlist_all == d};

        selectedFunction = ii;
        y = applyToRows(@(x)testfunction(x,'selectedfunction',selectedFunction),x);

        n = 14*ones(d,1);
        eftt = EFTT(@(x)testfunction(x,'selectedfunction',selectedFunction),d,'n',n,'tol',tol,'basis','legendre',...
            'TuckerRanks',maxTuckerranks,'ttranks',maxTTranks,'verb',0);
        EFTTerror(kk,ii) = norm(fevaluate(eftt,x) - y)/norm(y);
        EFTTeval(kk,ii) = eftt.numAllEvals;
        EFTTdof(kk,ii) = dof(eftt);
        EFTTdegrees{kk,ii} = degree(eftt);
        EFTTttranks{kk,ii} = ttrank(eftt);
        EFTTtuckerranks{kk,ii} = tuckerrank(eftt);

        fprintf('| %11s | %8.2e | %7.i | %6.i |   %3.i  |   %2.i   |     %2.i     |\n',funcname{ii},...
            EFTTerror(kk,ii),EFTTeval(kk,ii),EFTTdof(kk,ii),max(EFTTdegrees{kk,ii}),max(EFTTttranks{kk,ii}),...
            max(EFTTtuckerranks{kk,ii}))
    end
    if savedata
        save('data/20BenchFuncsAdapt','-regexp', '^(?!(x_all|dlist_all)$).')
    end
end


%error, eval, dof, degrees, ttranks, tuckerranks
tabledata = zeros(n_funcs,6);
for ii = 1:n_funcs
    tabledata(ii,1) = geomean(EFTTerror(:,ii));
    tabledata(ii,2) = mean(EFTTeval(:,ii));
    tabledata(ii,3) = mean(EFTTdof(:,ii));
    tabledata(ii,4) = max(mean(cell2mat(EFTTdegrees(:,ii)'),1));
    tempttrank = cell2mat(EFTTttranks(:,ii)');
    tabledata(ii,5) = max(mean(tempttrank(2:end-1,:),1));
    tabledata(ii,6) = max(mean(cell2mat(EFTTtuckerranks(:,ii)'),1));
end

if savedata
    save('data/20BenchFuncsAdapt','-regexp', saveExceptThem)
end
%load('data/20BenchFuncsAdapt')
fprintf('-----------Table: EFTT (vs FTcross) for 20 benchmark functions with adaptive n-----------\n')
fprintf('|   function  |   error  |   eval  |   dof  | degree | TTrank | Tuckerrank |\n')
for ii = 1:n_funcs
    fprintf('| %11s | %8.2e | %7.i | %6.i |   %3.i  |   %2.i   |     %2.i     |\n',funcname{ii},...
        tabledata(ii,1),round(tabledata(ii,2)),round(tabledata(ii,3)),round(tabledata(ii,4)),...
        round(tabledata(ii,5)), round(tabledata(ii,6)))
end