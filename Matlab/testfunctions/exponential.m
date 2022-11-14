function [y] = exponential(xx)
% x \in [-1, 1]
y = -exp(-sum(xx.^2)/2);
end
