function [y] = schaffer(xx)
% x \in [-100, 100]
xx = xx.*100;
temp = xx(2:end).^2+xx(1:end-1).^2;
y = sum(  0.5+(sin(sqrt(temp)).^2-0.5)./(1+0.001*temp).^2  );
end
