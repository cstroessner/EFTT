function F = Vals2LegCoeffsMat(n,m)
% maps function evaluations at m Chebyshev points to Legendre coefficients

if n < 2 ||  m < 2
    warning('n too small')
end
F = zeros(n,m);
[~,w] = chebpts(m);
cheb = @(i,n) cos((i-1).*pi/(n-1));
xx = cheb(1:m,m);

P = zeros(n,m);
P(1,:) = ones(1,m);
P(2,:) = xx;
for i = 3:n
    P(i,:) = (2*i-3)/(i-1) * xx.* P(i-1,:) - (i-2)/(i-1) * P(i-2,:);
end


%for i = 1:n
%    F(i,:) = (2*i-1)/2 * w.* P(i,:);
%end
F = P .* ((2*(1:n)' - 1)/2 * w);
end