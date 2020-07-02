function [P,P_agg,P_cluster,P_cluster_agg] = LR_Interaction(data,allgenes,cluster_label,L,R,TG_act,TG_inh)
% Input:
% data:             gene expression of size m x n (m genes and n cells)
% allgenes:         gene annotations
% cluster_label:    cell id (cluster id of all cells)
% L:                ligand (paired with R)
% R:                receptor (paired with L)
% TG_act:           target genes (activation)
% TG_inh:           target genes (inhibition)

% Output:
% P:             probability matrix of cell-cell interactions for all L-R pairs
% P_agg:         aggegration of cell-cell interations for all L-R pairs
% P_cluster:     cluster-cluster interaction matrix for all L-R pairs
% P_cluster_agg: aggegration of cluster-cluster interations for all L-R pairs 

No_pairs = length(L);
P = cell(No_pairs,1);
P_agg = zeros(length(cluster_label));

for i = 1:length(L)
    L_i = L{i};
    R_i = R{i};
    display([num2str(i) ':' '' L_i '-' R_i]);
    if nargin==5
        P1 = LR_Prob(data,allgenes,L_i,R_i);
        if max(P1(:)) > 0
            P1 = P1./max(P1(:));
        end
        P1(P1<=1e-6) = 0;
        P{i} = sparse(P1);
        P_agg = P_agg + P{i};
    elseif nargin==6
        %         TG_acti = TG_act(i);
        P1 = LRact_Prob(data,allgenes,L_i,R_i,TG_act);
        if max(P1(:)) > 0
            P1 = P1./max(P1(:));
        end

        P1(P1<=1e-6) = 0;
        P{i} = sparse(P1);
        P_agg = P_agg + P{i};
    elseif nargin==7
        P1 = LRboth_Prob(data,allgenes,L_i,R_i,TG_act,TG_inh);
        if max(P1(:)) > 0
            P1 = P1./max(P1(:));
        end

        P1(P1<=1e-6) = 0;
        P{i} = sparse(P1);
        P_agg = P_agg + P{i};
    end
end
if max(P_agg(:)) > 0
    P_agg = P_agg./max(P_agg(:));
end
P_agg = sparse(P_agg);

[P_cluster,P_cluster_agg] = P_agg_cluster(data,allgenes,P,L,R,cluster_label);
end



function [P,P_agg] = P_agg_cluster(data,allgenes,Pidv,L,R,cluster_label)
No_cluster = length(unique(cluster_label));
No_LR = size(Pidv,1);
P = cell(No_LR,1);
P_agg = zeros(No_cluster);
for i = 1:No_LR
    [~,~,zz1] = intersect(L(i),allgenes,'stable');
    [~,~,zz2] = intersect(R(i),allgenes,'stable');
    data1 = data(zz1,:);
    data2 = data(zz2,:);
    
    P0 = Pidv{i};
    P1 = zeros(No_cluster);
    for j = 1:No_cluster
        aa_idx = find(cluster_label==j);
        zz1 = length(find(data1(aa_idx)));
        for k = 1:No_cluster
            bb_idx = find(cluster_label==k);
            zz2 = length(find(data2(bb_idx)));
            cc = P0(cluster_label==j,cluster_label==k);
            % display(size(cc));
            % display([sum(cc(:)) zz1 zz2 zz1*zz2 log(zz1*zz2+1)]);
            if zz1*zz2 > 0
                P1(j,k) = sum(cc(:))./(zz1*zz2);
            end
        end
    end
    if max(P1(:)) > 0
        P1 = P1./max(P1(:));
    end

    P{i} = P1;
    P_agg = P_agg + P1;
end
if max(P_agg(:)) > 0
    P_agg = P_agg./max(P_agg(:));
end

end



function P = LR_Prob(data,allgenes,L,R)

[~,n] = size(data);

[~,L_idx,~] = intersect(allgenes,L,'stable');
[~,R_idx,~] = intersect(allgenes,R,'stable');

% Normalize L_data and R_data along gene
% Normalize L_data and R_data along cell
L_data = data(L_idx,:);
R_data = data(R_idx,:);

% display([allgenes{L_idx} '-' allgenes{R_idx}]);
display(['Expressed:' '' num2str(nnz(L_data)) '-' num2str(nnz(R_data))]);

P = zeros(n);

Lig_idx = find(L_data > eps);
Rec_idx = find(R_data > eps);

aa = L_data(Lig_idx);
bb = R_data(Rec_idx);

P1 = zeros(length(aa),length(bb));
for i = 1:length(aa)
    b = sum(exp(-1./( aa(i).*bb ) )  );
    
    for j = 1:length(bb)
        P1(i,j) = (exp(-1./(aa(i).*bb(j))))./b;
    end
end

P(Lig_idx,Rec_idx) = P1;
end


function P = LRact_Prob(data,allgenes,L,R,TG)

[~,n] = size(data);
[~,L_idx,~] = intersect(allgenes,L,'stable');
[~,R_idx,~] = intersect(allgenes,R,'stable');
[~,TG_idx,~] = intersect(allgenes,TG,'stable');
No_TGs = length(TG_idx);

L_data = data(L_idx,:);
R_data = data(R_idx,:);
TG_data = data(TG_idx,:);

if No_TGs > 1
    TG_all = mean(TG_data);
else
    TG_all = TG_data;
end

display(['Expressed:' '' num2str(nnz(L_data)) '-' num2str(nnz(R_data))]);

P = zeros(n);

Lig_idx = find(L_data);
Rec_idx = find(R_data);

aa = L_data(Lig_idx);
bb = R_data(Rec_idx);
cc = TG_all(Rec_idx); % gene expression on cells that have receptors

nn = length(aa);
nn1 = length(bb);
beta = exp(-1./cc);
Coef =zeros(nn,nn1);
for i = 1:nn
    alpha = (exp(-1./(aa(i).*bb)));
    Coef(i,:) = alpha./(alpha + beta);
end

P1 = zeros(length(aa),length(bb));
for i = 1:length(aa)
    alpha = (exp(-1./(aa(i).*bb)));
    b = sum(alpha.*Coef(i,:).*beta );
    
    for j = 1:length(bb)
        P1(i,j) = ((exp(-1./(aa(i).*bb(j)))).*Coef(i,j).*beta(j))./b;
    end
end

P(Lig_idx,Rec_idx) = P1;

end

function P = LRboth_Prob(data,allgenes,L,R,TG,TG1)
% Cell-cell interaction probability w.r.t L-R target activation (TG) and
% target inhibition (TG1)

[~,n] = size(data);
[~,L_idx,~] = intersect(allgenes,L,'stable');
[~,R_idx,~] = intersect(allgenes,R,'stable');
[~,TG_idx,~] = intersect(allgenes,TG,'stable');     % activator
[~,TG_idx1,~] = intersect(allgenes,TG1,'stable');   % inhibitor

No_TGs = length(TG_idx);

L_data = data(L_idx,:);
R_data = data(R_idx,:);
TG_data = data(TG_idx,:);
TG_data1 = data(TG_idx1,:);


TG_all = mean(TG_data);
TG_all1 = mean(TG_data1);

display(['Expressed:' '' num2str(nnz(L_data)) '-' num2str(nnz(R_data))]);


P = zeros(n);

Lig_idx = find(L_data>eps);
Rec_idx = find(R_data>eps);

aa = L_data(Lig_idx);
bb = R_data(Rec_idx);
cc = TG_all(Rec_idx); % gene expression on cells that have receptors
dd = TG_all1(Rec_idx); % gene expression on cells that have receptors

nn = length(aa);
nn1 = length(bb);
beta = exp(-1./cc);
beta1 = exp(-dd);

Coef =zeros(nn,nn1);
Coef1 =zeros(nn,nn1);

for i = 1:nn
    alpha = (exp(-1./(aa(i).*bb)));
    Coef(i,:) = alpha./(alpha + beta);
    Coef1(i,:) = alpha./(alpha + beta1);
end

P1 = zeros(length(aa),length(bb));
for i = 1:length(aa)
    alpha = (exp(-1./(aa(i).*bb)));
    b = sum(alpha.*Coef(i,:).*beta.*Coef1(i,:).*beta1);
    
    for j = 1:length(bb)
        P1(i,j) = ((exp(-1./(aa(i).*bb(j)))).*Coef(i,j).*beta(j).*Coef1(i,j).*beta1(j))./b;
    end
end

P(Lig_idx,Rec_idx) = P1;
end

