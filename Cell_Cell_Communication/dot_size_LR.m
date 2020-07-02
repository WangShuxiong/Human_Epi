function P = dot_size_LR(data,Lig,Rec,allgenes,cluster_labelr)
No_cluster = length(unique(cluster_labelr));
No_LR = length(Lig);

P = cell(No_cluster,1);

for i = 1:No_LR
    [~,~,L_idx] = intersect(Lig(i),allgenes);
    [~,~,R_idx] = intersect(Rec(i),allgenes);
    L_data = data(L_idx,:);
    R_data = data(R_idx,:);
    P0 = zeros(No_cluster);
    for j = 1:No_cluster
        for k = 1:No_cluster
            P0(j,k) = nnz(L_data(find(cluster_labelr==j)))*nnz(R_data(find(cluster_labelr==k)));
        end
    end
    P0 = log2(P0+1);
    P{i} = P0;
end


