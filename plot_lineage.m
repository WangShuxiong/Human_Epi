function plot_lineage(Lineage,No_cluster,cluster_label,Cell_dist,lgd,folder,mycolor,figname)

Nodesize = zeros(No_cluster,1);
for i = 1:No_cluster
    Nodesize(i) = length(find(cluster_label==i));
end
Nodesize = 30*Nodesize./max(Nodesize);

% plot cluster color on lineage tree
pred = Lineage;
rootedTree = digraph(pred(pred~=0),find(pred~=0));

plot(rootedTree,'Marker','o','MarkerSize',Nodesize,'NodeColor',mycolor(1:No_cluster,:),'NodeLabel',lgd);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
% title('Cluster on lineage')
box off;
axis off;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = [0 0 2 4];

% print([folder '\Lineage_Cluster_Color_' figname],'-depsc','-r300'); 
print([folder '\Lineage_cluster' figname],'-dpdf','-r300'); 



cmap = parula;
mymap = cmap(1:round(256*58./64),:);

ptimecolor = zeros(No_cluster,1);

for i = 1:No_cluster
     ptimecolor(i) = mean(Cell_dist(cluster_label==i));
%     ptimecolor(i) = mean(Ptime(find(cluster_label==i)));
end

pred = Lineage;
rootedTree = digraph(pred(pred~=0),find(pred~=0));
figure;
colormap(mymap);
plot(rootedTree,'Marker','o','MarkerSize',Nodesize,'NodeCData',ptimecolor, 'NodeColor','flat','NodeLabel',lgd);

% cb = colorbar;
% ax = gca;
% axpos = ax.Position;
% cpos = cb.Position;
% cpos(3) = 0.5*cpos(3);
% cb.Position = cpos;
% ax.Position = axpos;
% % lim = caxis
% % cb.Limits = lim;
% aa = cell(1,2);
% aa{1} = 'low';
% aa{2} = 'high';
% cb.TickLabels{1} = aa{1};
% cb.TickLabels{end} = aa{2};
% 
% for ii = 2:length(cb.TickLabels)-1
%     cb.TickLabels{ii} = [];
% end
% box off;

axis off;
% set(gca,'LineWidth',1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'FontName','Arial');
set(gca,'FontSize',12);


fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

fig.Units = 'Inches';
fig.Position = [0 0 2 4];
box off;

% title('Pseudotime on lineage')
print([folder '\Lineage_ptime_' figname],'-dpdf','-r300');  %,'-fillpage'
% print([folder '\Lineage_ptime_Color_' figname],'-depsc','-r300'); 
