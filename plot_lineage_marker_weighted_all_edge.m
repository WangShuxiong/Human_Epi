function plot_lineage_marker_weighted_all_edge(Lineage,No_cluster,cluster_label,Cell_dist,CC_adjacent,lgd,No_added_edges,figname,folder,data,marker,allgenes)

Nodesize = zeros(No_cluster,1);
for i = 1:No_cluster
    Nodesize(i) = length(find(cluster_label==i));
end
Nodesize = 30*Nodesize./max(Nodesize);


% plot cluster color on lineage tree
pred = Lineage;
% aa = pred(pred~=0);
% bb = find(pred~=0);
% Tw = zeros(length(aa),1);
% for i = 1:length(aa)
%     Tw(i) = WM(aa(i),bb(i));
% end
rootedTree = digraph(pred(pred~=0),find(pred~=0));
% Lwidth = 5*rootedTree.Edges.Weight./max(rootedTree.Edges.Weight);

% plot(rootedTree,'Marker','o','MarkerSize',Nodesize,'NodeColor',mycolor(1:No_cluster,:),...
%     'LineWidth',Lwidth,'ArrowSize',12,'EdgeColor',[0.69 0.77 0.87],...
%     'NodeLabel',lgd,'layout','force'); %,'WeightEffect','direct'

lt = plot(rootedTree); % lineage: lineage tree

xc = lt.XData;
yc = lt.YData;
% delete(lt);


% Add extra edges to the tree
adj_tree = zeros(length(lgd));
z = rootedTree.Edges.EndNodes;
for i = 1:length(z)
    adj_tree(z(i,1),z(i,2)) = 1;
    adj_tree(z(i,2),z(i,1)) = 1;
end

CC_adjacent1 = CC_adjacent;
CC_adjacent1(adj_tree == 0) = 0;

val = sort(CC_adjacent(adj_tree == 0));
val = val(val>0);

% No_added_edges = 3;
bb = zeros(size(CC_adjacent));
bb(adj_tree == 0) = CC_adjacent(adj_tree == 0);
bb(bb > val(No_added_edges)) = 0;

CC_adjacent1(adj_tree == 0) = bb(adj_tree == 0);

% figure;
bg = graph(CC_adjacent1,'upper');
Lwidth = 5*bg.Edges.Weight./max(bg.Edges.Weight);



[~,~,ia] = intersect(marker,allgenes,'stable');
No_mks = length(ia);
mk = allgenes(ia);

data_mk = data(ia,:);

aa1 = hot;
bb = bone;
zz = bb(176:end,:);
mycolormap = [zz;flip(aa1(96:2:end,:))];

for j = 1:No_mks
    mycolor = zeros(No_cluster,1);
    for i = 1:No_cluster
        mycolor(i) = mean(data_mk(j,cluster_label==i));
    end
    figure;
    colormap(mycolormap);
    
    aa = plot(bg,'Marker','o','MarkerSize',Nodesize,'NodeCData',mycolor,...
        'LineWidth',Lwidth,'EdgeColor',[0.69 0.77 0.87],...
        'NodeColor','flat','NodeLabel',[]); %,'WeightEffect','direct','layout','force'
    aa.XData = xc;
    aa.YData = yc;
    
    axis off;
    set(gca,'FontName','Arial');
    set(gca,'FontSize',12);
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    fig.Units = 'Inches';
    fig.Position = [0 0 3 4];
    
    % print([folder '\Lineage_Cluster_all_edges_' figname],'-depsc','-r300');
    print([folder '\' figname '_' mk{j}],'-dpdf','-r300');
end