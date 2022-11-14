function [y] = dixon(xx)
% x \in [-10, 10]
xx = 10*xx;
y = (xx(1)-1).^2+ sum((2*xx(2:end).^2 - xx(1:end-1)).^2.*(2:numel(xx)));
end
