function plot_pseudotime(latent,Cell_dist,folder,figname,method)
%% plot pseudotime on the latent space
cmap = parula;
mymap = cmap(1:round(60*256./64),:);
colormap(mymap);

scatter(latent(:,1),latent(:,2),10,Cell_dist,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);

% axis off;

cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;
% box on;
aa = cell(1,2);
aa{1} = '0';
aa{2} = '1';
cb.TickLabels{1} = aa{1};
cb.TickLabels{end} = aa{2};

for ii = 2:length(cb.TickLabels)-1
    cb.TickLabels{ii} = [];
end

xlabel([method num2str(1)]);
ylabel([method num2str(2)]);
set(gca,'FontName','Arial');
set(gca,'FontSize',12);

ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;

% Set the size of output fig
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

fig.Units = 'Inches';
fig.Position = [0 0 4 3];


% print([folder '\' figname '_' method],'-dpdf','-r300','-bestfit');
