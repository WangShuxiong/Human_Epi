function plot_mean_matrix(data,marker,allgenes,lgd,cluster_label,reorder,figname,folder)

No_cluster = length(unique(cluster_label));
[~,~,gene_idx] = intersect(marker,allgenes,'stable');

tg_mean = zeros(length(gene_idx),No_cluster);
data1 = data(gene_idx,:);

for i = 1:No_cluster
    tg_mean(:,i) = mean(data1(:,cluster_label==i),2);
end

tg_mean = tg_mean./max(tg_mean,[],2);

colormap(redbluecmap);
imagesc(tg_mean(:,reorder));
hold on;
box off;

% set(gca,'ColorScale','log')
cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;


xticks(1:length(lgd(reorder)));
xticklabels(lgd(reorder));
xtickangle(90);

yticks(1:length(marker));
yticklabels(marker);

set(gca,'FontName','Arial');
set(gca,'FontSize',12);

ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;


fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

fig.Units = 'Inches';
fig.Position = [0 0 6 10];

print([folder '\' figname],'-depsc','-r300');
