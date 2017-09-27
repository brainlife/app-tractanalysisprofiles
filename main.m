function [] = main()

switch getenv('ENV')
    case 'IUHPC'
        disp('loading paths for IUHPC')
        %addpath(genpath('/N/u/brlife/git/encode'))
        addpath(genpath('/N/u/brlife/git/vistasoft'))
        addpath(genpath('/N/u/brlife/git/jsonlab'))
        %addpath(genpath('/N/u/brlife/git/afq-master'))
    case 'VM'
        disp('loading paths for Jetstream VM')
        %addpath(genpath('/usr/local/encode'))
        addpath(genpath('/usr/local/vistasoft'))
        addpath(genpath('/usr/local/jsonlab'))
        %addpath(genpath('/usr/local/afq-master'))
end

% load my own config.json
config = loadjson('config.json');
%config.dt6 config.afq
%dt = dtiLoadDt6('/N/dc2/projects/lifebid/code/kitchell/app-tractanalysisprofiles/dtiinit/dti/dt6.mat');
dt = dtiLoadDt6(fullfile(config.dt6,'/dti/dt6.mat'));
load(config.afq);
%load('output.mat');
prop = config.diff_prop;

mkdir('images');

for ifg = 1:length(fg_classified)
    fg = fg_classified( ifg );

    % compute the core fiber from the fiber group (the tact profile is computed here)
    % [fa, md, rd, ad, cl, core] = dtiComputeDiffusionPropertiesAlongFG( fg, dt,[],[],200);
    [fa, md, rd, ad, cl, core] = dtiComputeDiffusionPropertiesAlongFG( fg, dt,[],[],100);
    
    % How to make a trct profile from a NIFTI file (such as from a run model)
    % nifti_file = niftiRead('path/to/nifti/file.nii.gz')
    % val = dtiComputeDiffusionPropertiesAlongFG( fg, nifti_file,[],[],200);
    
    % 3. Select a center portion fo the tract and show the FA and MD values
    % normally we only use for analyses the middle most reliable portion of the fiber.
    % nodesToPlot = 50:151;
    nodesToPlot = 25:76;
    
    h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'off');
    
    switch prop
        case 'fa'
            tract_profile = plot(fa(nodesToPlot),'color', [0.2 0.2 0.9],'linewidth',4)
            ylh = ylabel('Fractional Anisotropy');
            ylim = [0.00, 1.00];
            ytick = [0 .25 .5 .75];
        case 'md'
            tract_profile = plot(md(nodesToPlot),'color', [0.2 0.2 0.9],'linewidth',4)
            ylh = ylabel('Mean Diffusivity');
            ylim = [0.00, 1.00];
            ytick = [0 .25 .5 .75];
        case 'rd'
            tract_profile = plot(rd(nodesToPlot),'color', [0.2 0.2 0.9],'linewidth',4)
            ylh = ylabel('Radial Diffusivity');
            ylim = [0.00, 1.00];
            ytick = [0 .25 .5 .75];
        case 'ad'
            tract_profile = plot(ad(nodesToPlot),'color', [0.2 0.2 0.9],'linewidth',4)
            ylh = ylabel('Axial Diffusivity');
            ylim = [0.00, 2.00];
            ytick = [0 .5 1 1.5];
    end
    set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
        'xticklabel',{'Tract begin','Tract end'},'xlim',[0 50],'ylim',ylim,'Ytick',ytick,'Xtick',[0 50])
    Title_plot = title(fg.name);
    xlabel('Location on tract')
    
    
    saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_', prop)), 'png')
    saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_', prop)), 'eps')
    %clear()
    json.images(ifg).filename = strcat('images/',Title_plot.String,'_', prop, '.png');
    json.images(ifg).name = strcat(Title_plot.String);
    json.images(ifg).desc = strcat(Title_plot.String, ' tract analysis profile');
    
end
clf
savejson('', json, fullfile('images.json'));
end

