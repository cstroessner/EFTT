function tabledata = test_GenzFuncsAdapt(n_realizations,savedata)
% Test the Genz functions for EFTT with n choosen adaptively
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
                fprintf('|    function   |  method  |   error  |   eval  |   dof  | degree | TTrank | Tuckerrank |\n')
            end

            d_index = (dlist_all == d);
            w = w_all{d_index,jj}(kk,:);
            c = c_all{d_index,jj}(kk,:);
            x = x_all{d_index};
            y = applyToRows(@(x)testfunction(x,'w',w,'c',c,'selectedfunction',selectedFunction),x);

            n = 14*ones(d,1);
            eftt = EFTT(@(x)testfunction(x,'w',w,'c',c,'selectedfunction',selectedFunction),d,'n',n,'tol',tol,...
                'TuckerRanks',maxTuckerranks,'ttranks',maxTTranks,'verb',0,'basis','legendre');
            EFTTerror(kk,ii,jj) = norm(fevaluate(eftt,x) - y)/norm(y);
            EFTTeval(kk,ii,jj) = eftt.numAllEvals;
            EFTTdof(kk,ii,jj) = dof(eftt);
            EFTTdegrees{kk,ii,jj} = degree(eftt);
            EFTTttranks{kk,ii,jj} = ttrank(eftt);
            EFTTtuckerranks{kk,ii,jj} = tuckerrank(eftt);

            fprintf('| %13s |   EFTT   | %8.2e | %7.i | %6.i |   %3.i  |   %2.i   |     %2.i     |\n',funcname{jj},...
                EFTTerror(kk,ii,jj),EFTTeval(kk,ii,jj),EFTTdof(kk,ii,jj),max(EFTTdegrees{kk,ii,jj}),max(EFTTttranks{kk,ii,jj}),...
                max(EFTTtuckerranks{kk,ii,jj}))
        end
    end
    if savedata
        save('data/GenzFuncsAdapt','-regexp', saveExceptThem)
    end
end

%error, eval, dof, degrees, ttranks, tuckerranks
tabledata = zeros(n_funcs,length(dlist1),6);
for ii = 1:n_funcs
    for jj = 1:length(dlist1)
        tabledata(ii,jj,1) = geomean(EFTTerror(:,jj,ii));
        tabledata(ii,jj,2) = mean(EFTTeval(:,jj,ii));
        tabledata(ii,jj,3) = mean(EFTTdof(:,jj,ii));
        tabledata(ii,jj,4) = max(mean(cell2mat(EFTTdegrees(:,jj,ii)'),1));
        tempttrank = cell2mat(EFTTttranks(:,jj,ii)');
        tabledata(ii,jj,5) = max(mean(tempttrank(2:end-1,:),1));
        tabledata(ii,jj,6) = max(mean(cell2mat(EFTTtuckerranks(:,jj,ii)'),1));
    end
end

if savedata
    save('data/GenzFuncsAdapt','-regexp', saveExceptThem)
end
%load('data/GenzFuncsAdapt')
fprintf('-----------Table: EFTT (vs FTcross) for Genz functions with adaptive n-----------\n')
for ii = 1:n_funcs
    fprintf('|    function   |  d  |  method  |   error  |   eval  |   dof  | degree | TTrank | Tuckerrank |\n')
    for jj = 1:length(dlist1)
        d = dlist1(jj);
        if ii == 2
            d = dlist2(jj);
        end

        if jj == 4
            fprintf('| %13s | %3.i |   EFTT   | %8.2e | %7.i | %6.i |   %3.i  |   %2.i   |     %2.i     |\n',...
                funcname{ii},d,tabledata(ii,jj,1),round(tabledata(ii,jj,2)),round(tabledata(ii,jj,3)),round(tabledata(ii,jj,4)),round(tabledata(ii,jj,5)),...
                round(tabledata(ii,jj,6)))
        else
            fprintf('|               | %3.i |   EFTT   | %8.2e | %7.i | %6.i |   %3.i  |   %2.i   |     %2.i     |\n',...
                d,tabledata(ii,jj,1),round(tabledata(ii,jj,2)),round(tabledata(ii,jj,3)),round(tabledata(ii,jj,4)),round(tabledata(ii,jj,5)),...
                round(tabledata(ii,jj,6)))
        end
    end
    fprintf('------------------------------------------------------------------------\n')
end
