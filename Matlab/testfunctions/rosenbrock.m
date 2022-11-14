function [y] = rosenbrock(xx)
% x \in [-2.048, 2.048]
xx = xx.*2.048;
y = sum(100.*(xx(2:end)-xx(1:end-1).^2).^2+(1-xx(1:end-1)).^2);
end
