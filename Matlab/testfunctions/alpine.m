function [y] = alpine(xx)
% x \in [-10, 10]
xx = 10*xx;
y = sum(abs(xx.*sin(xx)+0.1.*xx));
end
