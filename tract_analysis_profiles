function tract_analysis_profiles(fg,metric,metricName,numnodes,ylab,ylim,ytick)
    h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'off');
    
    tract_profile = plot(metric,'color', [0.2 0.2 0.9],'linewidth',4)
    ylh = ylabel(sprintf('%s',ylab));
    set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
    'xticklabel',{'Tract begin','Tract end'},'xlim',[0 numnodes],'ylim',ylim,'Ytick',ytick,'Xtick',[0 numnodes])
    Title_plot = title(fg.name);
    xlabel('Location on tract')
    saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, sprintf('_%s',metricName))), 'png')
    saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, sprintf('_%s',metricName))), 'eps')
end
