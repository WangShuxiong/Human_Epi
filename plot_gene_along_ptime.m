function plot_gene_along_ptime(data,marker,Ptime,allgenes,folder)
[~,Ptime1] = sort(Ptime);
% plot gene expression along ptime
[mk,~,~] = intersect(allgenes,marker);

No_cell = size(data,2);
alpha = 0;
for i = 1:length(mk)
    [~,~,ia] = intersect(mk{i},allgenes,'stable');
    mk_val = data(ia,:);
    
    % fit curve for each branch
    polyorder1 = 5;
 
    x = 0:1/(No_cell-1):1;
    y = mk_val(Ptime1);
    x_idx = find(y>alpha);
    
    p_b1 = polyfit(x(x_idx),y(x_idx),polyorder1);
    f_b1 = polyval(p_b1,x(x_idx));
    
    % plot gene expression along ptime and branch
    c = linspace(0,1,No_cell);
%     figure(i);
    scatter(x,mk_val(Ptime1),10,c,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    
    hold on;
    plot(x(x_idx),f_b1,'--','LineWidth',2,'Color','k','MarkerSize',10);
    xticks([0,1]);   
    
    title(allgenes{ia});
    ax = gca;
    fig = gcf;
%     fig = figure(i);
    ax.XTickMode = 'manual';
    ax.YTickMode = 'manual';
    ax.ZTickMode = 'manual';
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    %
    % fig.PaperUnits = 'inches';
    % fig.PaperPosition = [0 0 4 2];
    %
    fig.Units = 'Inches';
    fig.Position = [0 0 4 3];
    
%     print([folder '\ptime_' allgenes{ia}],'-dpdf','-r300'); % -dtiff
    saveas(gcf,[folder '\ptime_' allgenes{ia}],'epsc');
end