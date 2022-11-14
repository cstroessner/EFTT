function fhat = fevaluate(ftt,x)

% x should be of size s * d, s is the sample size, d is dimension

% transform from EFTT to FTT
[s,d] = size(x);

n = degree(ftt);
r = ttrank(ftt);

batchsize = 1000;
nodes = [1:batchsize:s,s+1];
fhat = zeros(s,1);

% evaluation
%%% if you encounter a memory error here, decrease the batchsize
if ftt.basis == "chebyshev"
    cheb = @(i,x) cos((i)*acos(x));
    for k = 1:(length(nodes)-1)
        fhattemp = reshape(sum(core(ftt.Ahat,1).*reshape(cheb((0:(n(1)-1))',x(nodes(k):nodes(k+1)-1,1)'),1,n(1),1,[]),2),r(1),r(2),[]);
        for i=2:d
            temp = reshape(sum(core(ftt.Ahat,i).*reshape(cheb((0:(n(i)-1))',x(nodes(k):nodes(k+1)-1,i)'),1,n(i),1,[]),2),r(i),r(i+1),[]);
            fhattemp = pagemtimes(fhattemp,temp);
        end
        fhat(nodes(k):nodes(k+1)-1,1) = squeeze(fhattemp);
    end
elseif ftt.basis == "legendre"
    for k = 1:(length(nodes)-1)
        fhattemp=reshape(sum(core(ftt.Ahat,1).*reshape(myLegPoly(n(1),x(nodes(k):nodes(k+1)-1,1)),1,n(1),1,[]),2),r(1),r(2),[]);
        d = length(n);
        for i=2:d
            temp=reshape(sum(core(ftt.Ahat,i).*reshape(myLegPoly(n(i),x(nodes(k):nodes(k+1)-1,i)),1,n(i),1,[]),2),r(i),r(i+1),[]);
            fhattemp = pagemtimes(fhattemp,temp);
        end
        fhat(nodes(k):nodes(k+1)-1,1) = squeeze(fhattemp);
    end
end
end