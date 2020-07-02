function plot_lineage_weighted_all_edge(Lineage,No_cluster,cluster_label,Cell_dist,CC_adjacent,lgd,No_added_edges,figname,folder,mycolor)

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

aa = plot(bg,'Marker','o','MarkerSize',Nodesize,'NodeColor',mycolor(1:No_cluster,:),...
    'LineWidth',Lwidth,'EdgeColor',[0.69 0.77 0.87],...
    'NodeLabel',lgd,'layout','force'); %,'WeightEffect','direct'
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
print([folder '\Lineage1_' figname],'-dpdf','-r300'); 



cmap = parula;
mymap = cmap(1:58*256/64,:);
ptimecolor = zeros(No_cluster,1);

for i = 1:No_cluster
     ptimecolor(i) = mean(Cell_dist(cluster_label==i));
end
ptimecolor = ptimecolor - min(ptimecolor);
ptimecolor = ptimecolor./max(ptimecolor);

figure;
aa = plot(bg,'Marker','o','MarkerSize',Nodesize,'NodeCData',ptimecolor,...
    'NodeColor','flat','EdgeLabel',[],'LineWidth',Lwidth,'EdgeColor',[0.69 0.77 0.87],...
    'NodeLabel',lgd,'layout','force'); %,'WeightEffect','direct','EdgeLabel',Lwidth
aa.XData = xc;
aa.YData = yc;

% 
% cb = colorbar;
% ax = gca;
% axpos = ax.Position;
% cpos = cb.Position;
% cpos(3) = 0.5*cpos(3);
% cb.Position = cpos;
% ax.Position = axpos;
% 
% cb.Ticks = [0 1];
% %cb.TickDirection = 'out';
% % lim = caxis;
% % cb.Limits = 0.5*lim;
% aa = cell(1,2);
% aa{1} = '0';
% aa{2} = '1';
% cb.TickLabels{1} = aa{1};
% cb.TickLabels{end} = aa{2};
% 
% for ii = 2:length(cb.TickLabels)-1
%     cb.TickLabels{ii} = [];
% end



axis off;
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = [0 0 3 4];

% print([folder '\Lineage_Ptime_all_edges_' figname],'-depsc','-r300'); 
print([folder '\Lineage1_' figname],'-dpdf','-r300'); 
