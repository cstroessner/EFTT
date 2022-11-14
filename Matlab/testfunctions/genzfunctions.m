function [val] = genzfunctions(x,w,c)
global selectedFunction
d = size(x,2);
modified = 0;
if modified
    switch selectedFunction
        case 1
            func = @(x) cos(2*pi*w(1)+ sum(c.*x)); %oscillatory
        case 2
            func = @(x) prod(1./(c.^(-2)+ (x + w).^2));%product peak
        case 3
            func = @(x) (1+sum(c.*x)).^(-d-1);%corner peak
        case 4
            func = @(x) exp(-sum(c.^2.*(x-w).^2)); %Gaussian
        case 5
            func = @(x) exp(-sum(c.^2.*abs(x-w))); %C0
        case 6
            func = @(x) exp(sum(c.*x))*(x(1)<=w(1))*(x(2)<=w(2)); %discontinous
    end
else
    b = [9;7.25;1.85;7.03;20.4;4.3]';
    h = [1.5;2;2;1;2;2]';
    b = 10.^h.*b;
    
    c = c./norm(c,1)*b(selectedFunction)/(d^h(selectedFunction));
    switch selectedFunction
        case 1
            func = @(x) cos(2*pi*w(1)+ sum(c.*x)); %oscillatory
        case 2
            func = @(x) prod(1./(c.^(-2)+ (x + w).^2));%product peak
        case 3
            func = @(x) (1+sum(c.*x)).^(-d-1);%corner peak
        case 4
            func = @(x) exp(-sum(c.^2.*(x-w).^2)); %Gaussian
        case 5
            func = @(x) exp(-sum(c.^2.*abs(x-w))); %C0
        case 6
            func = @(x) exp(sum(c.*x))*(x(1)<=w(1))*(x(2)<=w(2)); %discontinous
    end
end
val = applyToRows(func,(x+1)./2);
global numfevals
numfevals = numfevals + numel(val);
end
