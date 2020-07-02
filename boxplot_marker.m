function boxplot_marker(data,allgenes,marker,cluster_labs,No_cluster,folder)
% Box plot for each gene along all clusters

% colormap jet;
% cmap1 = jet;
% mymap1 = cmap1(1:end,:);
% ncolor = size(mymap1,1);
% mycolor = mymap1(1:round(ncolor./No_cluster):1+ncolor,:);


zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];



group = cell(size(cluster_labs));
for i = 1:length(cluster_labs)
    group{i} = ['C' num2str(cluster_labs(i))];
end

ia = zeros(1,length(marker));
for i = 1:length(marker)
    for j = 1:length(allgenes)
        if strcmp(upper(marker{i}),upper(allgenes{j}))
            ia(i) = j;
        end
    end
end

% [~,ia,~] = intersect(allgenes,marker,'stable');
gname = allgenes(ia);

display(allgenes(ia));

MM = data(ia,:);

for i = 1:No_cluster
    cluster_notation{i} = ['C' num2str(i)];
end

n = length(ia);
for i = 1:n
    figure;
    boxplot(MM(i,:),group,'GroupOrder',cluster_notation,'Notch','on','PlotStyle','compact','Widths',0.9,'Colors',mycolor,...
        'LabelOrientation','horizontal'); 
    title(gname{i});
    print(['Results\' folder '\Markers\box_mk_' marker{i}],'-dpdf','-r300'); %'-dpdf',
end