% This script generates the data for table 1

initialize_QoI()

clear;
savedata = true;
mlist = [2,3,4];
tollist = [1e-3,1e-4];
maxTTranks = 25;
maxTuckerranks = 25;
load("testPoints.mat",'x_all','dlist_all')
load("cookiepde.mat",'fcookie2','fcookie3','fcookie4')

saveExceptThem = '^(?!(x_all|dlist_all|fcookie2|fcookie3|fcookie4|fcookie)$).';

fprintf('|  d  |   tol   |  method  |   error  |   eval  |   dof  |  time (min)  |  degree | TTrank | Tuckerrank |\n')
for ii = 1:length(mlist)
    m = mlist(ii);
    d = m^2;
    n = 14*ones(d,1);
    x = x_all{dlist_all == d};
    switch m
        case 2
            fcookie = @(x)fcookie2(x);
        case 3
            fcookie = @(x)fcookie3(x);
        case 4
            fcookie = @(x)fcookie4(x);
    end
    y = applyToRows(@(x)fcookie(x),x);

    for jj = 1:length(tollist)
        tol = tollist(jj);
        rng(0)
        tic
        eftt = EFTT(@(x)fcookie(x),d,'n',n,'tol',tol,'ifadapt',ones(d,1),...
            'TuckerRanks',maxTuckerranks,'ttranks',maxTTranks,'verb',0,'basis','legendre');
        EFTTtime(ii,jj) = toc;
        EFTTerror(ii,jj) = norm(fevaluate(eftt,x) - y)/norm(y);
        EFTTeval(ii,jj) = eftt.numAllEvals;
        EFTTdof(ii,jj) = dof(eftt);
        EFTTdegrees{ii,jj} = degree(eftt);
        EFTTttranks{ii,jj} = ttrank(eftt);
        EFTTtuckerranks{ii,jj} = tuckerrank(eftt);

        fprintf('| %3.i | %7.1e |   EFTT   | %8.2e | %7.i | %6.i | %10.2f   |   %3.i   |   %2.i   |     %2.i     |\n',d,tol,...
            EFTTerror(ii,jj),EFTTeval(ii,jj),EFTTdof(ii,jj),EFTTtime(ii,jj)/60,max(EFTTdegrees{ii,jj}),max(EFTTttranks{ii,jj}),max(EFTTtuckerranks{ii,jj}))

        if savedata
            save('data/cookiepde_adapt','-regexp', saveExceptThem)
        end
    end
end