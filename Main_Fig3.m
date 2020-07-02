% add folder to matlab path
addpath('Data');
addpath('CCS')
addpath('Results');
%% load data
load Data_Fig3.mat;

%% set folder to save all figures
folder = 'Results\Fig3';
folder1 = 'Results/Fig3';

%% Ligand receptor pairs in the list
[Lig,Rec] = ligand_receptor_in_list(allgenes,Lig,Rec);
target_up = intersect(allgenes,target_up);
target_down = intersect(allgenes,target_down);

%% Infer cell-cell (cluster-cluster) interation probability matrix via SoptSC
[P_cell_each_lr,P_cell_agg,P_cluster_each_lr,P_cluster_agg] = LR_Interaction(data,allgenes,cluster_lab,Lig,Rec,target_up,target_down);

%% Clustering L-R pairs
text_on = 1; % 0: not l-r names; 1: l-r names on;
NC_LR = 7;   % Number of clusters for l-r pairs
figname1 = 'Cluster_LR_Wnt';
cluster_lab_lr = cluster_LR(P_cluster_each_lr,NC_LR,Lig,Rec,figname1,folder,text_on);


%% plot average expression of selected genes in each cluster
figname1 = 'LR_mean';
gene_selected = [unique(Lig);unique(Rec);unique(target_up);unique(target_down)];
reorder = 1:length(lgd);
plot_mean_matrix(data,gene_selected,allgenes,lgd,cluster_lab,reorder,figname1,folder);


%% plot aggegration of cluster-cluster interations (aggegration of all l-r pair)
figname = 'Cluster-cluster-interaction-agg';
threshold = 0.3;
plot_sig_network(P_cluster_agg,cluster_lab,threshold,lgd,folder,figname,cluster_color);

%% dot size
P_dot_size = dot_size_LR(data,Lig,Rec,allgenes,cluster_lab);

%% plot dot-heatmap at cluster-level for all lig-rec pairs
figname = 'Dot_heatmap';
cluster_order_lr = 1:NC_LR;
reorder = 1:length(lgd);
plot_dot_heatmap(P_cluster_each_lr,P_dot_size,Lig,Rec,lgd,cluster_color,color_lr_all,reorder,cluster_lab_lr,cluster_order_lr,folder,figname);

%% plot cluster-cluster interactions based on the groups of ligand-receptor pairs
threshold = 0.1;
figname = 'Signaling_cluster';
cluster_order_lr = 1:NC_LR;
plot_sig_network_group_lr(P_cluster_each_lr,lgd,threshold,cluster_lab,cluster_lab_lr,cluster_order_lr,cluster_color,figname,folder);


%% plot cluster-cluster interations for a specific ligand-receptor pair
lr_number = 3; % plot the 3rd l-r pair

figname = [Lig{lr_number} '-' Rec{lr_number}];
threshold = 0.3;
plot_sig_network(P_cluster_each_lr{lr_number},cluster_lab,threshold,lgd,folder,figname,cluster_color);



