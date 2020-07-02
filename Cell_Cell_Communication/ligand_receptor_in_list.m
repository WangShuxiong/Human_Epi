function [Lig,Rec] = ligand_receptor_in_list(allgenes,lig,rec)

Lig = [];
Rec = [];
for i = 1:length(lig)
    [aa,~,~] = intersect(allgenes,[lig(i);rec(i)],'stable');
    if length(aa) == 2
        Lig = [Lig; lig(i)];
        Rec = [Rec; rec(i)];
    end
end


