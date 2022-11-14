function tabledata = test_SineIntegrationAdapt(n_repetitions,savedata)
% Test rank 2 sin function integration in Legendre basis
% The results in the paper are obtained by running the experiment 100 times
% and taking an average. You may simply set n_repetitions to 1 to get the
% figure quickly.

%savedata = false;
n_repetitions = 1;
selectedFunction = 'sine';
dlist = floor(linspace(10,450,10));
tollist = logspace(-13, -1, 10);
tollist = tollist([1,5,10]);


saveExceptThem = '^(?!(kk|ii|jj)$).';
for kk = 1:n_repetitions
    rng(kk)
    fprintf('-----------repetition = %i-----------\n',kk)
    fprintf('|  d  |   tol   |   error  |  eval |  dof  | degree |\n')
    for ii = 1:length(dlist)
        d = dlist(ii);
        I_true = imag(((exp(1i)-1)/1i)^d)*2^d;
        for jj = 1:length(tollist)
            tol = tollist(jj);

            eftt = EFTT(@(x)testfunction(x,'selectedfunction',selectedFunction),d,'n',6*ones(d,1),'tol',tol,'basis','legendre',...
                'TuckerRanks',2,'ttranks',2,'verb',0,'num_check',2);

            EFTTerror(kk,ii,jj) = abs(fintegrate(eftt) - I_true)/abs(I_true);
            EFTTeval(kk,ii,jj) = eftt.numAllEvals;
            EFTTdof(kk,ii,jj) = dof(eftt);
            EFTTdegrees{kk,ii,jj} = degree(eftt);
            EFTTttranks{kk,ii,jj} = ttrank(eftt);
            EFTTtuckerranks{kk,ii,jj} = tuckerrank(eftt);

            fprintf('| %3.i | %7.1e | %8.2e | %5.i | %5.i |   %2.i   |\n',d,tol,...
                EFTTerror(kk,ii,jj),EFTTeval(kk,ii,jj),EFTTdof(kk,ii,jj),max(EFTTdegrees{kk,ii,jj}));
        end
    end
    if savedata
        save('data/SineIntegrationAdapt','-regexp', saveExceptThem)
    end
end


%error, eval, dof, degrees
tabledata = zeros(length(dlist),length(tollist),4);
EFTTerror(EFTTerror==0) = min(EFTTerror(EFTTerror~=0),[],'all');
for ii = 1:length(dlist)
    for jj = 1:length(tollist)
        tabledata(ii,jj,1) = geomean(EFTTerror(:,ii,jj));
        tabledata(ii,jj,2) = mean(EFTTeval(:,ii,jj));
        tabledata(ii,jj,3) = mean(EFTTdof(:,ii,jj));
        tabledata(ii,jj,4) = max(mean(cell2mat(EFTTdegrees(:,ii,jj)'),1));
    end
end

if savedata
    save('data/SineIntegrationAdapt','-regexp', saveExceptThem)
end
% load('data/SineIntegrationAdapt')
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