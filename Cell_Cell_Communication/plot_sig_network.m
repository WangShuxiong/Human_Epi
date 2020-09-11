function plot_sig_network(P_cluster,cluster_label,threshold,lgd,folder,figname,mycolor)
% plot cell-cell and cluster-cluster signaling network based on the probability matrix
% Input:
%   -- Pidv: Cell-to-cell interaction probability of individual
%      ligand-receptor pair where Pidv{i} is the probability corresponds to
%      the ith ligand-receptor pair.
%   -- Pall: cell-to-cell interaction probability based on all
%      ligand-receptor pairs and their target genes.
%   -- threshold: restricting probabiltiy between cells less than threshold
%      to be zero.
%   -- folder: folder name where the results will be saved to.
%
%   Output:
%   -- cell-cell signaling network for individual ligand-receptor pair and
%   all ligand-receptor pairs.

No_cluster = length(unique(cluster_label));
zz = (0:No_cluster)./(No_cluster);
tickval = 0.5.*(zz(2:end) + zz(1:end-1));


adjacentM =P_cluster;
adjacentM = adjacentM./max(adjacentM(:));
adjacentM(adjacentM < threshold) = 0;

if max(adjacentM(:)) > 0
    bg = digraph(adjacentM);
    bg.Edges.LWidths = 5*bg.Edges.Weight/max(bg.Edges.Weight);
    
    % Set my edge color
    aa = bg.Edges.EndNodes;
    myedgecolor = zeros(size(aa,1),3);
    for i = 1:size(aa,1)
        for j = 1:No_cluster
            if aa(i,1) == j
                myedgecolor(i,:) = mycolor(j,:);
            end
        end
    end
    Gh = plot(bg,'Marker','o','MarkerSize',20,'Layout','circle','NodeLabel',lgd,...
        'NodeColor',mycolor,'EdgeColor',myedgecolor,'LineWidth',...
        bg.Edges.LWidths,'ArrowSize',10,'Layout','layered');
    
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    axis off;
    title(figname);
    
    colormap(mycolor);
    cb = colorbar('Ticks',tickval,'TickLabels',lgd,'FontSize',12);
    
    %     cb = colorbar;
    ax = gca;
    axpos = ax.Position;
    cpos = cb.Position;
    cpos(3) = 0.5*cpos(3);
    cb.Position = cpos;
    ax.Position = axpos;
    box off;
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    fig.Units = 'Inches';
    fig.Position = [0 0 6 4];
    print([folder '\' figname],'-dpdf','-r300');
end

