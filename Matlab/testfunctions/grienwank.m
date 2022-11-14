function [y] = grienwank(xx)
% x \in [-600, 600]
xx = 600*xx;
y = sum(xx.^2)/4000 - prod(cos(xx)./(1:numel(xx))) + 1;
end
