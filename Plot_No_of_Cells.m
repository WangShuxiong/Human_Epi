function No_cells = Plot_No_of_Cells(cluster_labs,lgd,folder,reorder,mycolor,figname)

No_cells1 = [];
No_cluster = length(unique(cluster_labs));
zz1 = unique(cluster_labs);
for i = 1:No_cluster
    No_cells1 = [No_cells1; length(find(cluster_labs==zz1(i)))];
end


No_cells = zeros(length(lgd),1);
No_cells(zz1) = No_cells1;


No_cells = No_cells(reorder)./sum(No_cells);
b = bar(No_cells);
b.FaceColor = 'flat';

ylim([0 1]);
% labels1 = string(b(1).YData);

zz = No_cells;

xtips1 = b.XData;
ytips1 = b.YData;
labels2 = string(lgd);
for i = 1:length(lgd)
    labels2(i) = [num2str(round(zz(i)*100)) '%'];
end

text(xtips1,ytips1,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')


for i = 1:length(lgd)
    b.CData(i,:) = mycolor(i,:);
end
grid on;
box off;
title('Percentage of Cells');
xticks(1:length(lgd));
xticklabels(lgd(reorder));
xtickangle(45);

ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

fig.Units = 'Inches';
fig.Position = [0 0 4 2.5];


% print([folder '\Percent_Cell_' figname],'-dpdf','-r300'); %'-dpdf',,'-fillpage'
print([folder '\Percent_Cell_' figname],'-depsc','-r300'); %'-dpdf',,'-fillpage'
