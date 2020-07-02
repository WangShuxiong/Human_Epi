function plot_count(data,cluster_label,No_cluster,lgd,folder)

% cmap1 = jet;
% mymap1 = cmap1(1:end,:);
% ncolor = size(mymap1,1);
% 
% mycolor = mymap1(1:round(ncolor./No_cluster):1+ncolor,:);
zzz = get(gca,'colororder');
mycolor = zeros(12,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
mycolor(12,:) = [0.5 1 0.7];



h = violinplot(data,cluster_label',1:No_cluster);
    for j = 1:length(h)
        h(j).ViolinColor = mycolor(j,:);
    end
    xlim([0 No_cluster+1]);
set(gca,'Xtick',1:No_cluster)
set(gca,'Xticklabel',lgd,'FontSize',12);
title('Count_per_Cell');
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% 
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% 
% ax.XTickMode = 'manual';
% ax.YTickMode = 'manual';
% ax.ZTickMode = 'manual';
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';


%end
%ax = gca;
%fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];

%ax.XTickMode = 'manual';
%ax.YTickMode = 'manual';
%ax.ZTickMode = 'manual';
%ax.XLimMode = 'manual';
%ax.YLimMode = 'manual';
%ax.ZLimMode = 'manual';
%fig.Units = 'Inches';
%fig.Position = [0 0 10 8];

print(['Results\' folder '\Count_per_Cell_latent'],'-dpdf','-r300'); %'-dpdf',

 