function fhat = DirectTTEval(Ahat,x,n)
% x should be of size m * d, m is the sample size, d is dimension

% transform evaluation tensor That to coeffient tensor Ahat
% d = That.d;
% Y=cell(d,1);
% for k=1:d
%     [sc1,sc2,sc3] = size(core(That,k));
%     Y{k} = reshape(Vals2ChebCoeffsMat(n(k))*reshape(permute(core(That,k),[2,1,3]),[sc2,sc1*sc3]),[n(k),sc1,sc3]);
%     if k==1
%         Y{k} = squeeze(Y{k});
%     end
% end
% Ahat = tt_tensor(Y);

% evaluation
[m,d] = size(x);
cheb = @(i,x) cos((i)*acos(x));
fhat = reshape(sum(core(Ahat,1).*reshape(cheb((0:(n(1)-1))',x(:,1)'),1,[],1,m),2),Ahat.r(1),[],m);
for i=2:d
    temp = reshape(sum(core(Ahat,i).*reshape(cheb((0:(n(i)-1))',x(:,i)'),1,[],1,m),2),Ahat.r(i),[],m);
    fhat = pagemtimes(fhat,temp);
end
fhat = squeeze(fhat);
end