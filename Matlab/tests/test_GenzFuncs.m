function [tabledata,tablecompare] = test_GenzFuncs(n_realizations,savedata)
% Test the Genz functions for DirectTT and EFTT with fixed n = 100
% The results in the paper are obtained by testing 30 different
% realizations and taking an average. You may simply set
% n_realizations to 1 to get the table quickly.

rng(0)
%savedata = false;
%n_realizations = 30;
tol = 1e-10;
maxTTranks = 50;
maxTuckerranks = 50;

n_funcs = 6;
funcname = {'Oscillatory','ProductPeak','CornerPeak','Gaussian','Continuous','Discontinuous'};

dlist1 = [20,50,100,200,300,400,500];
dlist2 = [5,10,15,20,25,30,50];

load("GenzRealizations.mat",'c_all','w_all')
load("testPoints.mat",'x_all','dlist_all')
saveExceptThem = '^(?!(c_all|w_all|x_all|dlist_all)$).';

for kk = 1:n_realizations
    fprintf('-----------realization = %i-----------\n',kk)
    for ii = 1:length(dlist1)
        for jj = 1:n_funcs
            selectedFunction = funcname{jj};
            d = dlist1(ii);
            if strcmpi(selectedFunction,'productpeak')
                d = dlist2(ii);
            end
            if jj==1
                fprintf('-----------d=%i for ProductPeak, d=%i for others-----------\n',dlist2(ii),dlist1(ii))
                fprintf('|    function   |  method  |   error  |   eval  |   dof  | TTrank | Tuckerrank | reduce eval(%%) | reduce storage(%%) |\n')
            end

            d_index = (dlist_all == d);
            w = w_all{d_index,jj}(kk,:);
            c = c_all{d_index,jj}(kk,:);
            x = x_all{d_index};
            y = applyToRows(@(x)testfunction(x,'w',w,'c',c,'selectedfunction',selectedFunction),x);

            n = 100*ones(d,1);
            eftt = EFTT(@(x)testfunction(x,'w',w,'c',c,'selectedfunction',selectedFunction),d,'n',n,'tol',tol,'ifadapt',zeros(d,1),...
                'TuckerRanks',maxTuckerranks,'ttranks',maxTTranks,'verb',0);
            EFTTerror(kk,ii,jj) = norm(fevaluate(eftt,x) - y)/norm(y);
            EFTTeval(kk,ii,jj) = eftt.numAllEvals;
            EFTTdof(kk,ii,jj) = dof(eftt);
            EFTTttranks{kk,ii,jj} = ttrank(eftt);
            EFTTtuckerranks{kk,ii,jj} = tuckerrank(eftt);


            ftt = FTT(@(x)testfunction(x,'w',w,'c',c,'selectedfunction',selectedFunction),d,'n',n,'tol',tol,'ttranks',maxTTranks,'verb',0);
            FTTerror(kk,ii,jj) = norm(fevaluate(ftt,x) - y)/norm(y);
            FTTeval(kk,ii,jj) = ftt.numAllEvals;
            FTTdof(kk,ii,jj) = dof(ftt);
            FTTttranks{kk,ii,jj} = ttrank(ftt);

            CompareEval(kk,ii,jj) = (1-EFTTeval(kk,ii,jj)/FTTeval(kk,ii,jj))*100;
            CompareDof(kk,ii,jj) = (1-EFTTdof(kk,ii,jj)/FTTdof(kk,ii,jj))*100;

            fprintf('| %13s |   EFTT   | %8.2e | %7.i | %6.i |   %2.i   |     %2.i     |     %6.1f     |       %6.1f      |\n',funcname{jj},...
                EFTTerror(kk,ii,jj),EFTTeval(kk,ii,jj),EFTTdof(kk,ii,jj),max(EFTTttranks{kk,ii,jj}),...
                max(EFTTtuckerranks{kk,ii,jj}),CompareEval(kk,ii,jj),CompareDof(kk,ii,jj))
            fprintf('|               | DirectTT | %8.2e | %7.i | %6.i |   %2.i   |            |                |                   |\n',...
                FTTerror(kk,ii,jj),FTTeval(kk,ii,jj),FTTdof(kk,ii,jj),max(FTTttranks{kk,ii,jj}))
        end
    end
    if savedata
        save('data/GenzFuncs','-regexp', saveExceptThem)
    end
end

%error, eval, dof, ttranks, tuckerranks
tabledata = zeros(n_funcs,length(dlist1),10);
for ii = 1:n_funcs
    for jj = 1:length(dlist1)
        tabledata(ii,jj,1) = geomean(EFTTerror(:,jj,ii));
        tabledata(ii,jj,3) = mean(EFTTeval(:,jj,ii));
        tabledata(ii,jj,5) = mean(EFTTdof(:,jj,ii));
        tempttrank = cell2mat(EFTTttranks(:,jj,ii)');
        tabledata(ii,jj,7) = max(mean(tempttrank(2:end-1,:),1));
        tabledata(ii,jj,9) = max(mean(cell2mat(EFTTtuckerranks(:,jj,ii)'),1));

        tabledata(ii,jj,2) = geomean(FTTerror(:,jj,ii));
        tabledata(ii,jj,4) = mean(FTTeval(:,jj,ii));
        tabledata(ii,jj,6) = mean(FTTdof(:,jj,ii));
        tempttrank = cell2mat(FTTttranks(:,jj,ii)');
        tabledata(ii,jj,8) = max(mean(tempttrank(2:end-1,:),1));
        tabledata(ii,jj,10) = 0;
    end
end

% reduce eval (%), reduce storage (%), increse accuracy (digit)
tablecompare = zeros(n_funcs,length(dlist1),3);
tablecompare(:,:,1) = 1 - tabledata(:,:,3)./tabledata(:,:,4);
tablecompare(:,:,2) = 1 - tabledata(:,:,5)./tabledata(:,:,6);
tablecompare(:,:,3) = log10(tabledata(:,:,2)./tabledata(:,:,1));

if savedata
    save('data/GenzFuncs','-regexp', saveExceptThem)
end
%load('data/GenzFuncs')
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
