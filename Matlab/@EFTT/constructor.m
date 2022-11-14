function eftt = constructor(eftt,testfunctionnD,d,varargin)
% Things may increase efficiency:
% 1. use catch warning in backslash rather than rcond
% no need 2. If the dof of tucker core is smaller than dof of TT decomposition, use
% Tucker core directly
% done 3. store as EFTT format, now it is TT format
% 4. TT decompostion of Tucker core
% done 5. create TT tensor using tt_tensor({cell of core}), now it is tt_random

% input
% fun: vectorized function handle x \in N*d \mapsto @(x)fun(x) \in N*1
% d: dimension of the function
% optional
% tol: tolerance of ACA and greedy2cross
% maxTuckerRanks:
% maxTTranks:
% n: vector of d * 1, number of parameters, default is [17;17;...;17]
% basis: "chebyshev": chebyshev interpolation on the chebyshev grid
%  or "legendre": legendre projection with coeffients computed by Clenshaw–Curtis quadrature
% ifadapt: vector of size d * 1, if adaptively choosing the size of grid, default is [1;1;...;1]
% ifchop: if computing the chopped degrees or not, default is [1;1;...;1]
% use_getTol: use getTol or not
% tol_aca: tolerance for ACA
% tol_chop: tolerance for chopping
% tol_ttcross: tolerance for greedy2cross
% num_check: number of coefficient to check for the Legendre basis
% projection

%% Initializations
if (~isempty(varargin))
    v1 = varargin{1};
    if (isa(v1, 'cell'))
        varargin=v1;
    end
end
vars = varargin;

% default values
basis = "chebyshev";
tol = 1e-10;
n = 17*ones(d,1); % size of initial Chebyshev grid
maxTuckerRanks = 20;
maxTTranks = 20;
ifadapt = ones(d,1);
use_getTol = 0;
tol_aca = tol;
tol_chop = tol;
tol_ttcross = tol;
num_check = 4;
verb = 1;


i=1;
while (i<length(vars))
    switch lower(vars{i})
        case 'tol'
            tol=vars{i+1};
            tol_aca = tol;
            tol_chop = tol;
            tol_ttcross = tol;
        case 'tuckerranks'
            maxTuckerRanks = vars{i+1};
        case 'ttranks'
            maxTTranks = max(vars{i+1},2);
        case 'ifadapt'
            ifadapt=vars{i+1};
        case 'basis'
            basis=vars{i+1};
        case 'n'
            n=vars{i+1};
        case 'use_gettol'
            use_getTol = vars{i+1};
        case 'tol_aca'
            tol_aca = vars{i+1};
        case 'tol_chop'
            tol_chop = vars{i+1};
        case 'tol_ttcross'
            tol_ttcross = vars{i+1};
        case 'num_check'
            num_check = vars{i+1};
        case 'verb'
            verb = vars{i+1};
        otherwise
            warning('Option %s was not recognized', vars{i});
    end
    i=i+2;
end


if basis == "legendre"
    n2m = @(n)ceil((n+1)./2); %@(n) n; 
    m2n = @(m)(2.*m-1); %@(m) m;
    m = n; %degree of Legendre
    n = m2n(m); %grid size
end

cn = zeros(d,1); % size of core tensor
U = cell(d,1); % approximation of span
TF = cell(d,1); % Tucker factor matrices
neval = 0;
% vectorize test function
testfunctionnD = @(x) applyToRows(testfunctionnD,x);
%% ACA
n = reshape(n,1,[]);
T = @(ind) testfunctionnD(cos((ind - 1) * pi./(n - 1)));

for i = 1:d
    %B = @(ind) zeros(size(ind,1),1);
    excepti = [1:i-1,i+1:d];
    
    mypivots = [];
    Ui = [];
    for ite = 1:maxTuckerRanks
        samplesize = min(ceil(geomean(n)/2),50);
        if ite == 1
            Tres = T;
        else
            Tres = @(ind)(T(ind) - B(ind));
        end
        maxvalue = 0;
        num_sampling = 0;
        
        % sample at most 5 times in the first iteration, sample once
        % afterwards
        while maxvalue == 0 && num_sampling <= 5
            if ite >= 2
                num_sampling = 5;
            end
            num_sampling = num_sampling + 1;
            randomIndex = generateRandomIndex(d,n,samplesize);
            M = Tres(randomIndex);
            if ifadapt(i) == 0
                neval = neval + samplesize + 2*samplesize*(ite-1);
            else
                neval = neval + samplesize + (ite-1)^2 + 2*samplesize*(ite-1);
            end
            [maxvalue, maxcolindex] = max(abs(M));
        end
        
        if verb > 0
        if basis == "chebyshev"
            fprintf('=ACA= d=%i, ite=%d, degree=%i, maxvalue=%3.3e, ACA#evals=%d\n', i, ite,n(i),maxvalue, neval);
        elseif basis == "legendre"
            fprintf('=ACA= d=%i, ite=%d, degree=%i, maxvalue=%3.3e, ACA#evals=%d\n', i, ite,n2m(n(i)),maxvalue, neval);
        end
        end
        
        if maxvalue < tol_aca && (ite>=2)
            break
        end
        
        maxval_list(ite,1) = maxvalue;
        [~,argmin_now] = min(maxval_list);
        
        if (ite>=15) && (ite - argmin_now >=5)
            break
        end
        
        if use_getTol == 1
            [reltol, tol] = getTol(M, eps, tol);
            tol_aca = tol;
            tol_chop = reltol;
            tol_ttcross = tol;
        end
        
        pp = randomIndex(maxcolindex,:);
        Ui = [Ui,MatrixUnfolding(T,i,(1:n(i))',pp(:,excepti))];
        neval = neval + n(i);
        
        mypivots = [mypivots;pp];
        
        if ifadapt(i)==1
            if basis == "chebyshev"
                
                while standardChop(Vals2ChebCoeffsMat(n(i)) * Ui(:,ite), tol_chop) == n(i) && n(i) <= 129
                    n(i) = 2*n(i) - 1;
                    T = @(ind) testfunctionnD(cos((ind - 1) * pi./(n - 1)));
                    tempUi = zeros(n(i)+1,ite);
                    for k = 1:ite
                        tempUi(:,k) = reshape([Ui(:,k),[MatrixUnfolding(T,i,(2:2:n(i))',mypivots(k,excepti));0]]',[],1);
                        neval = neval + (n(i) - 1)/2;
                    end
                    Ui = tempUi(1:n(i),:);
                    %Ui - MatrixUnfolding(T,i,(1:n(i))',mypivots(:,excepti))
                end
                
            elseif basis == "legendre"
                %while standardChop(Vals2LegCoeffsMat(m(i),n(i)) * Ui(:,ite), tol_chop) == m(i) && n(i) <= 5000
                while myLegChop(Vals2LegCoeffsMat(m(i),n(i)) * Ui(:,ite),'tol',tol_chop,'num_check',num_check) == m(i) && n(i) <= 105
                    n(i) = 2*n(i) - 1;
                    m(i) = n2m(n(i));
                    T = @(ind) testfunctionnD(cos((ind - 1) * pi./(n - 1)));
                    tempUi = zeros(n(i)+1,ite);
                    for k = 1:ite
                        tempUi(:,k) = reshape([Ui(:,k),[MatrixUnfolding(T,i,(2:2:n(i))',mypivots(k,excepti));0]]',[],1);
                        neval = neval + (n(i) - 1)/2;
                    end
                    Ui = tempUi(1:n(i),:);
                    %Ui - MatrixUnfolding(T,i,(1:n(i))',mypivots(:,excepti))
                end
            end
        end
        
        
        if maxvalue < tol_aca
            break
        end
        
        % compute cross matrix
        if ifadapt(i) == 0
            if ite == 1
                Aprev = T(mypivots);
                Anow = Aprev;
                neval = neval + 1;
            else
                Aprev = Anow;
                
                temp = mypivots(1:end-1,:);
                temp(:,excepti) = repmat(mypivots(end,excepti),ite-1,1);
                B1 = T(temp);
                
                temp = mypivots(1:end-1,:);
                temp(:,i) = mypivots(end,i);
                C1 = T(temp);
                D1 = T(pp);
                Anow = [Aprev, B1;C1',D1];
                neval = neval + 2*ite-1;
            end
        else
            Anow = MatrixUnfolding(T,i,mypivots(:,i),mypivots(:,excepti));
        end
        if rcond(Anow)<=eps
            break
        end
        B = @(ind) diag(MatrixUnfolding(T,i,ind(:,i),mypivots(:,excepti)) / Anow  * MatrixUnfolding(T,i,mypivots(:,i),ind(:,excepti)));
    end
    
    
    U{i} = Ui;
    cn(i) = size(Ui,2);
end
%% Implicit Tucker Approximation
maxcn = max(cn);
Im = zeros(maxcn,d);
for i = 1:d
    [Q,~] = qr(U{i},0);
    I = DEIM(Q);
    TF{i} = Q/Q(I,:);
    Im(1:cn(i),i) = I;
end
%% TT approximation of Tucker core
if maxcn==1
    corevalue = testfunctionnD(cos((Im - 1) * pi./(n - 1)));
    chat = tt_tensor(num2cell([ones(d-1,1);corevalue]));
    neval = neval + 1;
    evalttapprox = 1;
    swp = 1;
else
    [chat,~,swp,~,~,~,~,evalttapprox]=greedy2_cross(cn,@(ind) testfunctionnD( cos((Im(bsxfun(@plus,ind,maxcn*(0:d-1)))-1)*pi./(n-1)) ),tol_ttcross,'nswp',maxTTranks-1,'vec',true,'verb',verb);
    neval = neval + evalttapprox;
end
%% return values


if basis == "chebyshev"
    for i = 1:d
        TF{i} = Vals2ChebCoeffsMat(n(i))*TF{i};
    end
elseif basis == "legendre"
    for i = 1:d
        TF{i} = Vals2LegCoeffsMat(m(i),n(i))*TF{i};
    end
end

eftt.chat = chat;
eftt.TF = TF;

eftt.basis = basis;
eftt.numAllEvals = neval;
eftt.numTTEvals = evalttapprox;
eftt.numTTswp = swp;
end