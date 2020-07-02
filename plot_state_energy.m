function plot_state_energy(E_plus_state,order,lgd,mycolor,figname,folder)
%% bar plot
No_clusterb = length(order);

x = 1:No_clusterb;
fig = figure;
ax = axes('parent', fig);
p = plot(ax, x, E_plus_state(order));
p(1).LineWidth = 2;
p(1).Marker = 'o';
p(1).LineStyle = '--';
p(1).MarkerSize = 10;
p(1).Color = 'k';
% p(1).MarkerFaceColor = 'w';
p(1).MarkerEdgeColor = 'k';

xticks(1:No_clusterb);
xticklabels(lgd(order));
ylabel(ax, 'Energy of States')

hold on;
for ik = 1:No_clusterb
    scatter(x(ik),E_plus_state(order(ik)),60,mycolor(ik,:),'filled');
    hold on;
end
% xlabel(ax, 'Number of input data points')
% title(ax,'Comparison data for solving A\\b on different numbers of workers');

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

