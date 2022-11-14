addpath(pwd)
cd TT-Toolbox
addpath(pwd)
cd cross
addpath(pwd)
cd ..
cd core
addpath(pwd)
cd ..
cd ..
cd chebfun
addpath(pwd)
cd ..

clear
rng(0)
dlist_all = [3,4,5,6,7,8,9,10,15,16,20,25,30,40,50,60,70,80,90,100,200,300,400,500]';
for ii = 1:length(dlist_all)
    d = dlist_all(ii);
    x_all{ii,1} = 2*rand(10000,d)-1;
end

b = [9;7.25;1.85;7.03;20.4;4.3]';
h = [1.5;2;2;1;2;2]';
b = 10.^h.*b;
n_realizations = 30;
n_funcs = 6;
for ii = 1:length(dlist_all)
    for jj = 1:n_funcs
        d = dlist_all(ii);
        w_all{ii,jj} = rand(n_realizations,d);
        c_all{ii,jj} = rand(n_realizations,d);
        c_all{ii,jj} = c_all{ii,jj}./sum(c_all{ii,jj},2).*b(jj)./(d^h(jj));
    end
end

save('tests/testPoints','dlist_all','x_all')
save('tests/GenzRealizations','dlist_all','w_all','c_all')

addpath('tests/')
addpath('testfunctions/')
addpath('figure')
addpath('data')