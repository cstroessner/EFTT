function randomIndex = generateRandomIndex(d,n,num)

randomIndex = zeros(num,d);
for myindex = 1:d
randomIndex(:,myindex) = randi(n(myindex),num,1);
end