function plot_KRT1014(lgd,lgd_cass,order,data_mk1_mean_C,data_mk2_mean_C,folder,figname)
No_clusterb = length(lgd);
% Population Dynamics
x = 1:No_clusterb;
fig = figure;
ax = axes('parent', fig);

p = plot(ax, x, data_mk1_mean_C(order), ...
     x,data_mk2_mean_C(order));
 
% p = plot(ax, x, data_mk1_mean_C(order));

p(1).LineWidth = 2;
p(1).Marker = 'o';
p(1).LineStyle = '-';
p(1).MarkerSize = 10;
% p(1).MarkerFaceColor = 'w';
% p(1).MarkerEdgeColor = 'W';

p(2).LineWidth = 2;
p(2).LineStyle = ':';
p(2).Marker = 'v';
p(2).MarkerSize = 10;
% p(2).MarkerFaceColor = 'w';
% p(2).MarkerEdgeColor = 'W';

xticks(1:No_clusterb);
xticklabels(lgd(order));
ylabel(ax, 'Average Expression')
% xlabel(ax, 'Number of input data points')
% title(ax,'Comparison data for solving A\\b on different numbers of workers');
legend(lgd_cass,'Location','best');

% hold on;
% % for i = 1:No_clusterb
% scatter(x,data_mk1_mean_C(order),60,mycolor(order,:),'filled','o','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
% hold on;
% scatter(x,data_mk2_mean_C(order),60,mycolor(order,:),'filled','v','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
% end
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;

box off;
grid on;


fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = [0 0 5 3];

print([folder '\' figname],'-dpdf','-r300'); %'-dpdf',,'-fillpage'



