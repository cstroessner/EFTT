function ftt = constructor(ftt,testfunctionnD,d,varargin)
% input
% fun: vectorized function handle x \in N*d \mapsto @(x)fun(x) \in N*1
% d: dimension of the function

if (~isempty(varargin))
    v1 = varargin{1};
    if (isa(v1, 'cell'))
        varargin=v1;
    end
end
vars = varargin;

% default values
tol = 1e-10;
n = 100*ones(1,d); % size of initial Chebyshev grid
maxTTranks = 20;
verb = 0;


i=1;
while (i<length(vars))
    switch lower(vars{i})
        case 'n'
            n=reshape(vars{i+1},1,[]);
        case 'tol'
            tol=vars{i+1};
        case 'ttranks'
            maxTTranks = max(vars{i+1},2);
        case 'verb'
            verb = vars{i+1};
        otherwise
            warning('Option %s was not recognized', vars{i});
    end
    i=i+2;
end

% vectorize test function
testfunctionnD = @(x) applyToRows(testfunctionnD,x);

[That,~,ftt.numTTswp,~,~,~,~,ftt.numAllEvals]=greedy2_cross(n,@(ind) testfunctionnD(cos((ind - 1) * pi./(n - 1))),tol,'nswp',maxTTranks-1,'vec',true,'verb',verb);
%TTranks = That.r';
%dof = numel(That.core);

Y=cell(d,1);
for k=1:d
    [sc1,sc2,sc3] = size(core(That,k));
    Y{k} = reshape(Vals2ChebCoeffsMat(n(k))*reshape(permute(core(That,k),[2,1,3]),sc2,[]),[],sc1,sc3);
    if k==1
        Y{k} = squeeze(Y{k});
    end
end
ftt.Ahat = tt_tensor(Y);
end