function plot_lineage_weighted(Lineage,No_cluster,cluster_label,Cell_dist,CC_adjacent,lgd,figname,folder,mycolor)

Nodesize = zeros(No_cluster,1);
for i = 1:No_cluster
    Nodesize(i) = length(find(cluster_label==i));
end
Nodesize = 30*Nodesize./max(Nodesize);


% plot cluster color on lineage tree
WM = CC_adjacent + CC_adjacent';
pred = Lineage;
aa = pred(pred~=0);
bb = find(pred~=0);
Tw = zeros(length(aa),1);
for i = 1:length(aa)
    Tw(i) = WM(aa(i),bb(i));
end
rootedTree = digraph(pred(pred~=0),find(pred~=0),Tw);
Lwidth = 5*rootedTree.Edges.Weight./max(rootedTree.Edges.Weight);

% figure;
% plot(rootedTree,'Marker','o','MarkerSize',20,'NodeColor',mycolor(1:No_cluster,:),'NodeLabel',[]);
plot(rootedTree,'Marker','o','MarkerSize',Nodesize,'NodeColor',mycolor(1:No_cluster,:),...
    'LineWidth',Lwidth,'ArrowSize',12,'EdgeColor',[0.69 0.77 0.87],...
    'NodeLabel',lgd,'layout','force'); %,'WeightEffect','direct'

% box off;
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% title('Lineage Tree')
axis off;
set(gca,'FontName','Arial');
set(gca,'FontSize',12);

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = [0 0 3 4];


% print([folder '\Lineage_Cluster_Weighted_' figname],'-depsc','-r300'); 
print([folder '\Lineage_' figname],'-dpdf','-r300'); 


cmap = parula;
mymap = cmap(1:58*256/64,:);

ptimecolor = zeros(No_cluster,1);

for i = 1:No_cluster
     ptimecolor(i) = mean(Cell_dist(cluster_label==i));
%     ptimecolor(i) = mean(Ptime(find(cluster_label==i)));
end
ptimecolor = ptimecolor - min(ptimecolor);
ptimecolor = ptimecolor./max(ptimecolor);

pred = Lineage;
rootedTree = digraph(pred(pred~=0),find(pred~=0),Tw);
figure;
colormap(mymap);
plot(rootedTree,'Marker','o','MarkerSize',Nodesize,'NodeCData',ptimecolor,...
    'NodeColor','flat','NodeLabel',[],...
    'LineWidth',Lwidth,'ArrowSize',12,'EdgeColor',[0.69 0.77 0.87],...
    'NodeLabel',lgd,'layout','force');  % ,'WeightEffect','direct'

% cb = colorbar;
% % cb.Location = 'eastoutside';
% ax = gca;
% axpos = ax.Position;
% cpos = cb.Position;
% cpos(3) = 0.5*cpos(3);
% cb.Position = cpos;
% ax.Position = axpos;
% 
% cb.Ticks = [0 1];
% aa = cell(1,2);
% aa{1} = '0';
% aa{2} = '1';
% cb.TickLabels{1} = aa{1};
% cb.TickLabels{end} = aa{2};
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

% print([folder '\Lineage_Ptime_Weighted_' figname],'-depsc','-r300'); 
print([folder '\Lineage_' figname '1'],'-dpdf','-r300'); 
