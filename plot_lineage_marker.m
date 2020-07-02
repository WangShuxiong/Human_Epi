function plot_lineage_marker(data,allgenes,marker,Lineage,cluster_label,lgd,folder,figname)
No_cluster = length(unique(cluster_label));
[~,~,ia] = intersect(marker,allgenes,'stable');
No_mks = length(ia);
mk = allgenes(ia);

data_mk = data(ia,:);

pred = Lineage;
rootedTree = digraph(pred(pred~=0),find(pred~=0));
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
    plot(rootedTree,'Marker','o','MarkerSize',15,'NodeCData',mycolor,...
        'NodeColor','flat','NodeLabel',[]); % ,'NodeLabel',lgd
    axis off;
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    fig.Units = 'Inches';
    fig.Position = [0 0 2 4];
    
%     print([folder '\' figname '_' mk{j}],'-depsc','-r300');
    print([folder '\' figname '_' mk{j}],'-dpdf','-r300');
end