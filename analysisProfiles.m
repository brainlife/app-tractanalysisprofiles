function [imgnum,json] = analysisProfiles(data,fg,name,ytitle,ylim,ytick,numnodes,imgnum)

h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'off');
tract_profile = plot(data,'color', [0.2 0.2 0.9],'linewidth',4)
ylh = ylabel(ytitle);
ylim = ylim;
ytick = ytick;
set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
    'xticklabel',{'Tract begin','Tract end'},'xlim',[0 numnodes],'ylim',ylim,'Ytick',ytick,'Xtick',[0 numnodes])
Title_plot = title(fg.name);
xlabel('Location on tract')
saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_', name)), 'png')
saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_', name)), 'eps')
imgnum = imgnum+1;
json.images(imgnum).filename = strcat('images/',Title_plot.String,'_', name,'.png');
json.images(imgnum).name = strcat(Title_plot.String);
json.images(imgnum).desc = strcat(Title_plot.String, ' tract analysis profile');
clf;