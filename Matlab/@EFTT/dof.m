function dofs = dof(eftt)
dofs = numel(eftt.chat.core);
for ii = 1:eftt.chat.d
    dofs = dofs + numel(eftt.TF{ii});
end
end