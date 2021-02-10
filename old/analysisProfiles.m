function [] = analysisProfiles(data,fg,name,ytitle,ylim,ytick,numnodes,units)

h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'off');
tract_profile = plot(data,'color', [0.2 0.2 0.9],'linewidth',4)
ylh = ylabel({ytitle;strcat('(',units,')')});
ylim = ylim;
ytick = ytick;
set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
    'xticklabel',{'Tract RAS','Tract LPI'},'xlim',[0 numnodes],'ylim',ylim,'Ytick',ytick,'Xtick',[0 numnodes])
Title_plot = title(strrep(fg.name,'.','_'),'Interpreter','None');
Title_name = Title_plot.String;
xlabel('Location on tract')
saveas(tract_profile, fullfile('images/', strcat(strrep(Title_plot.String,' ','_'), '_', name)), 'png')
saveas(tract_profile, fullfile('images/', strcat(strrep(Title_plot.String,' ','_'), '_', name)), 'eps')
clf;
