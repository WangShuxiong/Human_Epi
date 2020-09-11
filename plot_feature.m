function plot_feature(data,mk,allgenes,latent,folder)
%This function plot gene expression along cell subpopulations
%   Input:
%  --   gene_set: genes to be plot specified by users
%  --   all_genes: all genes from the single cell data set
%  --   data: a m*n single cell data matrix with m rows(genes) and n columns(cells)
%  --   latent: low dimension space induced from cell-cell transition matrix
%
%   Output:
%           figures showing gene expression among cell subpopulations
%           on 2D-projection.
%

% MM(MM > 1.9) = 1.9;
mk = intersect(mk,allgenes,'stable');
gene_set_no = length(mk);

%% Marker genes expression on each subpopulation
aa1 = hot;
bb = bone;
zz = bb(176:end,:);
mymap = [zz;flip(aa1(96:2:end,:))];

for ik = 1:gene_set_no
    [~,~,gene_set_idx] = intersect(mk{ik},allgenes,'stable');
    MM = data(gene_set_idx,:);

    figure;
    colormap(mymap);
    %          colormap redbluecmap;
    scatter(latent(:,1),latent(:,2),10,MM,'filled','MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.8);
    
    box off;
    %     box off;
    %     set(gca,'LineWidth',1.5);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    axis off;
    set(gca,'FontName','Arial');
    set(gca,'FontSize',12);
    
    title(mk{ik});
    %         ax = gca;
    %         cb = colorbar;
    %         ax = gca;
    %         axpos = ax.Position;
    %         cpos = cb.Position;
    %         cpos(3) = 0.5*cpos(3);
    %         cb.Position = cpos;
    %         ax.Position = axpos;
    %         cb.TickLabels{1} = 0;
    %         cb.TickLabels{end} = 1;
    %
    %     for ii = 2:length(cb.TickLabels)-1
    %         cb.TickLabels{ii} = [];
    %     end
    
    ax = gca;
    cb = colorbar;
    axpos = ax.Position;
    cpos = cb.Position;
    cpos(3) = 0.5*cpos(3);
    cb.Position = cpos;
    ax.Position = axpos;
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    fig.Units = 'Inches';
    fig.Position = [0 0 4 3];
    box off;
    print([folder '\Feature_' mk{ik}],'-depsc','-r300'); %'-dpdf','-fillpage'

%      print([folder '\Feature_' mk{ik}],'-dpdf','-r300','-bestfit'); %'-dpdf',,'-bestfit'
end
