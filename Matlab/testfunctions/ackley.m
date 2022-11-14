function [y] = ackley(xx)
% x \in [−32.768, 32.768]
xx = 32.768*xx;
A = 20;B=0.2;C=2*pi;
y = -A*exp(-B*sqrt(mean(xx.^2))) - exp(mean(cos(C*xx)))+A+exp(1);
end
