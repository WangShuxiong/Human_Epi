function plot_eigengap(eigenvalues,folder)

scatter(1:min([30 size(eigenvalues,1)]),eigenvalues(1:min([30 size(eigenvalues,1)])),20,'filled');
grid on;
box off;
% set(gca,'LineWidth',1.5);
% xlabel('i');
% ylabel('Eigenvalue of Graph Laplacian');
yticks(0:0.2:1);
xticks(1:3:30)
set(gca,'FontName','Arial');
set(gca,'FontSize',12);

ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 4 3];

fig.Units = 'Inches';
fig.Position = [0 0 4 2.5];


print([folder '\EigenGap'],'-dpdf','-r300');
end
