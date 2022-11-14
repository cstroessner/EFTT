function [y] = qing(xx)
% x \in [0, 500]
xx = (xx+1)*250;
y = sum((xx.^2 - (1:numel(xx))).^2);
end
