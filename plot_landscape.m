function plot_landscape(latent_umap,Energy,cluster_labs,mycolor,folder,figname,method)
z = Energy;
x = latent_umap(:,1);
y = latent_umap(:,2);

xlin = linspace(min(x),max(x),100);
ylin = linspace(min(y),max(y),100);
[X,Y] = meshgrid(xlin,ylin);
% [X,Y] = meshgrid(x,y);
% f = scatteredInterpolant(x,y,z);
% Z = f(X,Y);

Z = griddata(x,y,z,X,Y,'cubic'); %'linear', 'natural','cubic'

%
% figure
% colormap redbluecmap;
% colormap jet;
mesh(X,Y,Z) %interpolated
% surfc(X,Y,Z);

% surf(X,Y,Z);
axis tight; 
% shading interp;

cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;


%% plot cells with cluster labs
hold on;

No_cluster = length(unique(cluster_labs));
for i = 1:No_cluster
    scatter3(x(cluster_labs==i),y(cluster_labs==i), z(cluster_labs==i),...
        20,mycolor(i,:),'filled','MarkerEdgeAlpha',0.9,'MarkerFaceAlpha',0.9);
% p = plot3(x(cluster_labs == i),y(cluster_labs == i),z(cluster_labs == i),'.','MarkerSize',5); %nonuniform
% p.Color = mycolor(i,:);

hold on;
end
% legend(lgd,'FontSize',12,'Location','best');%,'Orientation','horizontal');
% set(gca,'FontName','Arial');
set(gca,'xticklabels',[]);
set(gca,'yticklabels',[]);

box off;
xlabel([method num2str(1)]);
ylabel([method num2str(2)]);

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
fig.Position = [0 0 6 4];
box off;

print([folder '\' figname '_' method],'-dpdf','-r300'); %'-dpdf','-fillpage'


