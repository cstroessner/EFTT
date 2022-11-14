function fhat = fevaluate(eftt,x)

% x should be of size s * d, s is the sample size, d is dimension

% transform from EFTT to FTT
[s,d] = size(x);
Y=cell(d,1);

for k=1:d
    [sc1,sc2,sc3] = size(core(eftt.chat,k));
    Y{k} = reshape(eftt.TF{k}*reshape(permute(core(eftt.chat,k),[2,1,3]),sc2,[]),[],sc1,sc3);
    if k==1
        Y{k} = squeeze(Y{k});
    end
end

Ahat=tt_tensor(Y);

n = degree(eftt);
batchsize = 1000;
nodes = [1:batchsize:s,s+1];
fhat = zeros(s,1);

% evaluation
%%% if you encounter a memory error here, decrease the batchsize
if eftt.basis == "chebyshev"
    cheb = @(i,x) cos((i)*acos(x));
    for k = 1:(length(nodes)-1)
        fhattemp = reshape(sum(core(Ahat,1).*reshape(cheb((0:(n(1)-1))',x(nodes(k):nodes(k+1)-1,1)'),1,n(1),1,[]),2),Ahat.r(1),Ahat.r(2),[]);
        for i=2:d
            temp = reshape(sum(core(Ahat,i).*reshape(cheb((0:(n(i)-1))',x(nodes(k):nodes(k+1)-1,i)'),1,n(i),1,[]),2),Ahat.r(i),Ahat.r(i+1),[]);
            fhattemp = pagemtimes(fhattemp,temp);
        end
        fhat(nodes(k):nodes(k+1)-1,1) = squeeze(fhattemp);
    end
elseif eftt.basis == "legendre"
    for k = 1:(length(nodes)-1)
        fhattemp=reshape(sum(core(Ahat,1).*reshape(myLegPoly(n(1),x(nodes(k):nodes(k+1)-1,1)),1,n(1),1,[]),2),Ahat.r(1),Ahat.r(2),[]);
        d = length(n);
        for i=2:d
            temp=reshape(sum(core(Ahat,i).*reshape(myLegPoly(n(i),x(nodes(k):nodes(k+1)-1,i)),1,n(i),1,[]),2),Ahat.r(i),Ahat.r(i+1),[]);
            fhattemp = pagemtimes(fhattemp,temp);
        end
        fhat(nodes(k):nodes(k+1)-1,1) = squeeze(fhattemp);
    end
end
end