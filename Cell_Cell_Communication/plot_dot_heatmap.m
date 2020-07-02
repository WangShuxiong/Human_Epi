function [CC_single_M,CC_single_M1,xnames,ynames] = plot_dot_heatmap(P0,P1,Lig,Rec,lgd,mycolor,color_all,reorder,cluster_lab_lr,cluster_order_lr,folder,figname)

LR_order = [];
NC_LR = length(unique(cluster_lab_lr));
for i = 1:NC_LR
    LR_order = [LR_order;find(cluster_lab_lr==cluster_order_lr(i))];
end
cluster_lab_lr1 = cluster_lab_lr(LR_order);
P0 = P0(LR_order);
P1 = P1(LR_order);
Lig = Lig(LR_order);
Rec = Rec(LR_order);

% cluster_color for l-r pair
mycolor_lr = color_all(1:round(256./NC_LR):end,:);


No_LR = length(Lig);
No_cluster = length(lgd);

aa = [1:No_cluster^2]';
bb = repmat(aa,No_LR,1);
x = bb(:);

y = zeros(size(x));
for i = 1:No_LR
    y((i-1)*No_cluster^2 + 1:i*No_cluster^2) = No_LR + 1 - i;
end

CC_single_M = zeros(No_LR,No_cluster^2);    % dot color
CC_single_M1 = zeros(No_LR,No_cluster^2);   % dot size
for j = 1:No_LR
    zz0 = P0{j}';
    zz0 = zz0./max(zz0(:));
    CC_single_M(j,:) = reshape(zz0(reorder,reorder),1,No_cluster^2);
    zz = P1{j}';
    CC_single_M1(j,:) = reshape(zz(reorder,reorder),1,No_cluster^2);
end

% dot size
zz = CC_single_M1';
dot_size = zz(:);
% scale dot size
dot_size = 60.*dot_size./max(dot_size);

% dot color scale
zz = CC_single_M';
dot_color = zz(:);

colormap(redbluecmap);
scatter(x,y,dot_size +1,dot_color,'filled');
xlim([0 No_cluster^2+1]);
ylim([0 No_LR+1]);


ynames = cell(No_LR,1);
for i = 1:No_LR
    ynames{i} = [Lig{No_LR-i+1} '-' Rec{No_LR-i+1}];
end

xnames1 = repmat(lgd(reorder),No_cluster,1);
xnames0 = cell(size(xnames1));
lgd0 = lgd(reorder);
for i = 1:length(lgd)
    for j = 1:length(lgd)
        xnames0{j+(i-1)*length(lgd)} = lgd0{i};
    end
end

xnames = cell(size(xnames1));
for i = 1:length(xnames)
    xnames{i} = [xnames0{i} ' --> ' xnames1{i}];
end

yticks(1:length(ynames));
yticklabels(ynames);

xticks(1:length(xnames));
xticklabels(xnames1);
xtickangle(90);

cm = zeros(No_cluster^2,3);
for i = 1:No_cluster
    if i == 1
        cm(1:No_cluster,:) = ones(No_cluster,1).*mycolor(i,:);
    else
        cm(No_cluster*(i-1)+1:No_cluster*i,:) = ones(No_cluster,1).*mycolor(i,:);
    end
end
ax = gca;
for i = 1:No_cluster^2
    ax.XTickLabel{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        cm(i,:), ax.XTickLabel{i});
end



cm1 = zeros(No_LR,3);
for i = 1:No_LR
    for j = 1:NC_LR
        if cluster_lab_lr1(i) == cluster_order_lr(j)
            cm1(No_LR-i+1,:) = mycolor_lr(j,:);
        end
    end
end

% cm1 = cm1(LR_order,:);

for i = 1:No_LR
    ax.YTickLabel{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        cm1(i,:), ax.YTickLabel{i});
end


ax.TickDir = 'out';
ax.LineWidth = 1;
% box off;
% box on;
grid on;

% cb = colorbar('Ticks',tickval,'TickLabels',lgd,'FontSize',12);
cb = colorbar;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.3*cpos(3);
cb.Position = cpos;
ax.Position = axpos;

% Set the size of output fig
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = [0 0 10 8]; 
% print([folder '\Dot_heatmap__' figname],'-dpdf','-r300');%,'-fillpage','-dpsc'
% print([folder '\Dot_heatmap__' figname],'-dsvg','-r300');%,'-fillpage','-dpsc','-dsvg'

print([folder '\' figname],'-depsc','-r300');

