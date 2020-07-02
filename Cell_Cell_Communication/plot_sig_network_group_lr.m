function plot_sig_network_group_lr(Pidv_cluster,lgd,threshold,cluster_label,cluster_lab_lr,cluster_order_lr,mycolor,figname,folder)
NC_LR = length(unique(cluster_lab_lr));
No_cluster = length(lgd);

%% plot cell-cell signaling network based on cluster labels
zz = (0:No_cluster)./(No_cluster);
tickval = 0.5.*(zz(2:end) + zz(1:end-1));

Nodesize = zeros(No_cluster,1);
for i = 1:No_cluster
    Nodesize(i) = length(find(cluster_label==i));
end
% Nodesize = log2(100*Nodesize+1);
Nodesize = 20.*Nodesize./max(Nodesize) + 3;

for i = 1:NC_LR
    zz_idx = find(cluster_lab_lr == cluster_order_lr(i));
    
    P_cluster = zeros(No_cluster);
    for j = 1:length(zz_idx)
        zz1 =  Pidv_cluster{zz_idx(j)};
        if max(zz1(:)) > 0
            P_cluster = P_cluster +zz1./max(zz1(:));
        end
    end
    P_cluster = P_cluster./max(P_cluster(:));
    P_cluster(P_cluster < threshold) = 0;
    
    adjacentM = P_cluster;
    if max(adjacentM(:)) > 0
        %         adjacentM = adjacentM./max(adjacentM(:));
        %         bg1 = biograph(adjacentM);
        %         view(bg1);
        
        bg = digraph(adjacentM);
        bg.Edges.LWidths = 5*bg.Edges.Weight/max(bg.Edges.Weight);
        % Set my edge color
        aa = bg.Edges.EndNodes;
        myedgecolor = zeros(size(aa,1),3);
        for i1 = 1:size(aa,1)
            for jj = 1:No_cluster
                if aa(i1,1) == jj
                    myedgecolor(i1,:) = mycolor(jj,:);
                end
            end
        end
        
        figure;
        Gh = plot(bg,'Marker','o','MarkerSize',Nodesize,'NodeLabel',[],...
            'NodeColor',mycolor,'EdgeColor',myedgecolor,'LineWidth',...
            bg.Edges.LWidths,'ArrowSize',10,'Layout','circle'); % [0.690196 0.768627 0.870588]'force' 'Layout','layered'
        % 'Layout','circle', ,'XData',XY(:,1),'YData',XY(:,2)
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        box off;
        axis off;
        title(['L-R Group-' num2str(i)],'fontsize',12);
        
        colormap(mycolor);
        cb = colorbar('Ticks',tickval,'TickLabels',lgd,'FontSize',12);
        
        %         cb = colorbar;
        ax = gca;
        axpos = ax.Position;
        cpos = cb.Position;
        cpos(3) = 0.5*cpos(3);
        cb.Position = cpos;
        ax.Position = axpos;
        
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig.Units = 'Inches';
        fig.Position = [0 0 4 3];
        
%         print([folder '\' figname '_group_' num2str(i)],'-dpdf','-r300','-bestfit');
        print([folder '\' figname '_group_' num2str(i)],'-depsc','-r300');
    end
end
