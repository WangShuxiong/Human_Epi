function gene_idxv = GC_htmp_DE_lgd(data,allgenes,cluster_labs,topn,lgd,figname,folder,folder1)
%% new data visualization
No_clusterr = length(unique(cluster_labs));
gene_idxv = [];
cluster_order = [];

[No_gene,No_cell] = size(data);
% calculat mean of gene expression
gene_mean = zeros(No_gene,No_clusterr);
gene_DE_score = zeros(No_gene,No_clusterr);

% gene_value_idx = zeros(No_gene,1);
for i = 1:No_clusterr
    gene_mean(:,i) = mean(data(:,find(cluster_labs==i))');
    cluster_order = [cluster_order;find(cluster_labs==i)];
end
[~,gene_value_idx] = max(gene_mean,[],2);

% compute DE-score for each gene
for i = 1:No_clusterr
    zz = abs(gene_mean(:,i).*ones(1,No_clusterr) - gene_mean);
%     gene_DE_score(:,i) = mean(zz,2); % mean
    gene_DE_score(:,i) = sum(zz,2); % sum
end


% topn markers for each cluster based on DE score
gclusters = [];
gscore = [];
for i = 1:No_clusterr
    zz_idx = find(gene_value_idx == i);
    zz_DEscore = gene_DE_score(zz_idx,i);
    [zzvalue,zz1] = sort(zz_DEscore,'descend');
    topn1 = min([topn length(zz_idx)]);
    gene_idxv = [gene_idxv; zz_idx(zz1(1:topn1))];
    gclusters = [gclusters;i.*ones(topn1,1)];
    gscore = [gscore;zzvalue(1:topn1)];
end



GL500 = [gene_idxv gclusters gscore];
T = array2table(GL500,'RowNames',allgenes(gene_idxv),'VariableNames',{'Gene_indx','Cluster','DE_Score'});
writetable(T,[folder1 '/DE_Genes' num2str(topn) '_' figname '.csv'],'WriteRowNames',true,'WriteVariableNames',true,'Delimiter',',');  


% for i = 1:No_clusterr
%     zz_DEscore = gene_DE_score(:,i);
%     [~,zz1] = sort(zz_DEscore,'descend');
%     gene_idxv = [gene_idxv; zz1(1:topn)];
% end


datav = data(gene_idxv,cluster_order);
% datav = data_reduce1(Gene_labels_topnr(OGI,1),CGI);
% colormap redbluecmap;
% imagesc(datav); 
% figure;
idata = datav;
kk = 2;
center = mean(idata,kk);
scale = std(idata, 0,kk);
tscale = scale;
%=Check for zeros and set them to 1 so not to scale them.
scale(tscale == 0) = 1;
%== Center and scale the data
idata = bsxfun(@minus, idata, center);
sdata = bsxfun(@rdivide, idata, scale);
thresh = 3;
colormap redbluecmap;
clims = [-thresh thresh];
imagesc(sdata,clims);
set(gca,'xtick',[]);
set(gca,'ytick',[]);




if size(datav,1) < 200
yticks(1:size(datav,1));
yticklabels(allgenes(gene_idxv));
end


No_cells_inC = [];
for i = 1:No_clusterr
    No_cells_inC = [No_cells_inC; length(find(cluster_labs==i))];
end
xtkval = cumsum(No_cells_inC);
xtkval1 = zeros(size(xtkval));
for i = 1:No_clusterr
    if i==1
        xtkval1(i) = 0.5.*No_cells_inC(i);
    else
        xtkval1(i) = xtkval(i-1) + 0.5.*No_cells_inC(i);
    end
end

hold on;
% draw lines between clusters
yy = length(gene_idxv);
for i = 1:No_clusterr-1
     line([xtkval(i) xtkval(i)],[0 yy+0.5],'color','w',...
        'LineStyle','-','LineWidth',2);
    hold on;
end

xticks(xtkval1);
xticklabels(lgd);
xtickangle(45);
cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 4 3];
fig.Units = 'Inches';
fig.Position = [0 0 6 10];
box off;

print([folder '\HeatMap_Top' num2str(topn) '_' figname],'-depsc','-r300'); %'-dpdf', ,'-fillpage'

