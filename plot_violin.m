function plot_violin(data,allgenes,Marker,cluster_label,No_cluster,lgd,cluster_order,folder,mycolor)

cluster_labs1 = zeros(size(cluster_label));
for i = 1:No_cluster
    cluster_labs1(find(cluster_label==cluster_order(i))) = i;
end

for ik = 1:length(Marker)
    [~,ia,~] = intersect(allgenes,Marker{ik},'stable');
    figure(ik);
    h = violinplot(data(ia,:),cluster_labs1,1:No_cluster);
    for j = 1:length(h)
        h(j).ViolinColor = mycolor(j,:);
    end
    ylabel(Marker{ik});
    % hAxes = gca;
    % hAxes.XRuler.Axle.LineStyle = 'none';
    
    grid on;
    set(gca,'xtick',[]);
    axis off;
    
    set(gca,'ytick',[]);
    set(get(gca,'YLabel'),'visible','on')
    
    set(gca,'Xtick',1:No_cluster)
    set(gca,'Xticklabel',lgd(cluster_order));
    
    ax = gca;
    fig = gcf;
    ax.XTickMode = 'manual';
    ax.YTickMode = 'manual';
    ax.ZTickMode = 'manual';
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    
    fig.Units = 'Inches';
    fig.Position = [0 0 6 1];
    
%     print([folder '\violin_' Marker{ik}],'-dpdf','-r300'); %'-dpdf','-depsc','-r300'
    print([folder '\violin_' Marker{ik}],'-depsc','-r300'); %'-dpdf','-depsc','-r300'

end
