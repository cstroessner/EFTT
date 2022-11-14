function [y] = rastrigin(xx)
% x \in [-5.12, 5.12]
xx = xx.*5.12;
A = 10;
y = A*numel(xx) + sum(xx.^2 - A.*cos(2*pi.*xx));
end
