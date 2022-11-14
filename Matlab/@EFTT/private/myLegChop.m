function cutoff = myLegChop(coef,varargin)

if (~isempty(varargin))
    v1 = varargin{1};
    if (isa(v1, 'cell'))
        varargin=v1;
    end
end
vars = varargin;

% default values
tol = 1e-10;
num_check  = 4;

i=1;
while (i<length(vars))
    switch lower(vars{i})
        case 'tol'
            tol=vars{i+1};
        case 'num_check'
            num_check=vars{i+1};
        otherwise
            warning('Option %s was not recognized', vars{i});
    end
    i=i+2;
end

coef = abs(coef);
temp = prod( coef(end-num_check+1:end) < tol);
if temp == 0
    cutoff = length(coef);
elseif temp == 1
    cutoff = find(coef >= tol);
    if isempty(cutoff)
        cutoff = num_check;
    else
    cutoff = cutoff(end) + num_check;
    end
end
end