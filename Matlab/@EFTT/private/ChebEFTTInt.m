function I_est = ChebEFTTInt(EFTT)
%chat = EFTT{1};
%TF = EFTT{2};
d = EFTT{1}.d;
Y=cell(d,1);
for k=1:d
    [sc1,sc2,sc3] = size(core(EFTT{1},k));
    Y{k} = reshape(EFTT{2}{k}*reshape(permute(core(EFTT{1},k),[2,1,3]),sc2,sc1*sc3),[],sc1,sc3);
    if k==1
        Y{k} = squeeze(Y{k});
    end
end

Ahat=tt_tensor(Y);
I_est = core(Ahat,1);
I_est = reshape(squeeze(I_est(:,1,:)),1,[]);
for i = 2:d
    temp = core(Ahat,i);
    temp = squeeze(temp(:,1,:));
    I_est = reshape(I_est*temp,1,[]);
end
I_est = 2^d * I_est;
end