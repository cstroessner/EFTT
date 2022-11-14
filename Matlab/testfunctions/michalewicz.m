function [y] = michalewicz(xx)
% x \in [0, pi]
xx = (xx+1)*pi/2;
m = 10;
y = -sum(sin(xx).*sin(xx.^2.*(1:numel(xx))/pi).^(2*m));
end
