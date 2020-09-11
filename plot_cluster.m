function plot_cluster(latent,cluster_label,No_cluster,lgd,folder,reorder,method,figname,mycolor)


for ik = 1:No_cluster
    scatter(latent(find(cluster_label==reorder(ik)),1),latent(find(cluster_label==reorder(ik)),2),10,mycolor(ik,:),...
        'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
end
box off;
xlabel([method num2str(1)]);
ylabel([method num2str(2)]);
legend(lgd(reorder),'FontSize',12,'Location','eastoutside');%,'Orientation','horizontal');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;

% axis off;
% axis off;
% Set the size of output fig
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = [0 0 5 3];
box off;

print([folder '\' figname '_' method],'-depsc','-r300'); %'-dpdf','-fillpage'
% print([folder '\' figname '_' method],'-dpdf','-r300'); %'-dpdf','-fillpage'
