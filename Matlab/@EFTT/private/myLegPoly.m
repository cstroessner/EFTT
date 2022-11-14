function P=myLegPoly(n,x)
m = size(x,1);
P = zeros(n,m);
P(1,:) = 1;
P(2,:) = x;
for i = 3:n
    P(i,:) = (2*i-3)/(i-1) * x' .* P(i-1,:) - (i-2)/(i-1) * P(i-2,:);
end
end