function [val] = testfunction(x,varargin)

i=1;
while (i<length(varargin))
    switch lower(varargin{i})
        case 'w'
            w=varargin{i+1};
            d = length(w);
        case 'c'
            c = varargin{i+1};
        case 'selectedfunction'
            selectedFunction = varargin{i+1};
        otherwise
            warning('Option %s was not recognized', varargin{i});
    end
    i=i+2;
end

switch lower(selectedFunction)
    case 1
        val = ackley(x);
    case 2
        val = alpine(x);
    case 3
        val = dixon(x);
    case 4
        val = exponential(x);
    case 5
        val = grienwank(x);
    case 6
        val = michalewicz(x);
    case 7
        val = piston(x);
    case 8
        val = qing(x);
    case 9
        val = rastrigin(x);
    case 10
        val = rosenbrock(x);
    case 11
        val = schaffer(x);
    case 12
        val = schwefel(x);
    case 13
        val = borehole(x);
    case 14
        val = otlcircuit(x);
    case 15
        val = robot(x);
    case 16
        val = wingweight(x);
    case 17
        val = fried(x);
    case 18
        val = grlee09(x);
    case 19
        val = detpep108d(x);
    case 20
        val = detpep10exp(x);
    case 'sine'
        val = sin(sum(x+1)/2);
    case 'oscillatory'
        x = (x+1)/2;
        val = cos(2*pi*w(1)+ sum(c.*x));
    case 'productpeak'
        x = (x+1)/2;
        val = prod(1./(c.^(-2)+ (x + w).^2));
    case 'cornerpeak'
        x = (x+1)/2;
        val = (1+sum(c.*x)).^(-d-1);
    case 'gaussian'
        x = (x+1)/2;
        val = exp(-sum(c.^2.*(x-w).^2));
    case 'continuous'
        x = (x+1)/2;
        val = exp(-sum(c.^2.*abs(x-w)));
    case 'discontinuous'
        x = (x+1)/2;
        val = exp(sum(c.*x))*(x(1)<=w(1))*(x(2)<=w(2));
    otherwise
        val = 1;
end

global numfevals
numfevals = numfevals + numel(val);
end
