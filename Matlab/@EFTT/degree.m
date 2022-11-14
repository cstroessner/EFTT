function n = degree(eftt)
n = zeros(eftt.chat.d,1);
for ii = 1:eftt.chat.d
   n(ii) = size(eftt.TF{ii},1);
end
end