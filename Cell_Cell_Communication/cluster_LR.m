function cluster_lab_lr = cluster_LR(Pidv_cluster,NC_LR,Lig,Rec,figname,folder,text_on)

zz1 = parula;
No_cluster = size(Pidv_cluster{1},1);
No_LR = length(Lig);
M = zeros(No_LR,No_cluster^2);

for i = 1:No_LR
    aa = Pidv_cluster{i}';
    %      M(i,:) = aa(:)./max(aa(:));
    M(i,:) = aa(:);
end

% set up the initial
% flag = 2;
% [W0,H0] = nndsvd(M,NC_LR,flag);
%
% [W,H] = nnmf(M,NC_LR,'algorithm','mult','w0',W0,'h0',H0);
% [~,cluster_lab_lr] = max(W,[],2);

% % latent_pca = pca(M','NumComponents',2);
% % latent_pca = pca(W','NumComponents',2);

InitY = pca(M','NumComponents',2);
if No_LR > 100
    perplexity = 20;
else
    perplexity = 6;
end

if NC_LR <= perplexity
    perplexity = NC_LR;
end
latent_pca = tsne(M,'Algorithm','exact','InitialY',InitY,'Perplexity',perplexity);


Z = linkage(latent_pca,'average');
cluster_lab_lr = cluster(Z,'maxclust',NC_LR,'Depth',3);

mycolor_lr = zz1(1:round(256./NC_LR):end,:);

for ik = 1:NC_LR
    scatter(latent_pca(cluster_lab_lr==ik,1),latent_pca(cluster_lab_lr==ik,2),50,mycolor_lr(ik,:),'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
end
lgd = cell(NC_LR,1);
for i = 1:NC_LR
    lgd{i} = ['C' num2str(i)];
end
legend(lgd,'FontSize',12,'Location','eastoutside');%,'Orientation','horizontal');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;

set(gca,'xtick',[]);
set(gca,'ytick',[]);

xlabel('tSNE1');
ylabel('tSNE2');


lgd1 = {};
for i = 1:length(Lig)
    lgd1{i} = [Lig{i} '\_' Rec{i}];
end

if text_on == 1
    text(latent_pca(:,1),latent_pca(:,2),lgd1,'FontSize',9);
end

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

fig.Units = 'Inches';
fig.Position = [0 0 6 4];
box off;

% print([folder '\' figname],'-depsc','-r300'); %'-dpdf','-fillpage'
print([folder '\' figname],'-dpdf','-r300','-fillpage'); %'-dpdf','-fillpage'

