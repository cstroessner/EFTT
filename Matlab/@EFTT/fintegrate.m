function I_est = fintegrate(eftt)


d = eftt.chat.d;
Y=cell(d,1);
for k=1:d
    [sc1,sc2,sc3] = size(core(eftt.chat,k));
    Y{k} = reshape(eftt.TF{k}*reshape(permute(core(eftt.chat,k),[2,1,3]),sc2,sc1*sc3),[],sc1,sc3);
    if k==1
        Y{k} = squeeze(Y{k});
    end
end
Ahat=tt_tensor(Y);

if eftt.basis == "legendre"
    
    I_est = core(Ahat,1);
    I_est = reshape(squeeze(I_est(:,1,:)),1,[]);
    for i = 2:d
        temp = core(Ahat,i);
        temp = squeeze(temp(:,1,:));
        I_est = reshape(I_est*temp,1,[]);
    end
    I_est = 2^d * I_est;
    
elseif eftt.basis == "chebyshev"
    n = degree(eftt);
    core1 = core(Ahat,1);
    I_est = 2*reshape(squeeze(core1(:,1,:)),1,[]);
    for ik = 3:n(1)
        I_est = I_est + reshape(squeeze(core1(:,ik,:)),1,[])*((-1)^(ik-1) + 1)/(1-(ik-1)^2);
    end
    for i = 2:d
        corei = core(Ahat,i);
        temp = 2*squeeze(corei(:,1,:));
        for ik = 3:n(i)
            temp = temp + squeeze(corei(:,ik,:))*((-1)^(ik-1) + 1)/(1-(ik-1)^2);
        end
        I_est = reshape(I_est*temp,1,[]);
    end
end

end