function [] = main()

switch getenv('ENV')
    case 'IUHPC'
        disp('loading paths for IUHPC')
        %addpath(genpath('/N/u/hayashis/BigRed2/git/encode'))
        addpath(genpath('/N/u/hayashis/BigRed2/git/vistasoft'))
        addpath(genpath('/N/u/hayashis/BigRed2/git/jsonlab'))
        %addpath(genpath('/N/u/hayashis/BigRed2/git/afq-master'))
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

dt = dtiLoadDt6(fullfile(config.dt6,'/dti/dt6.mat'));
load(config.afq);

for ifg = 1:length(fg_classified)
    fg = fg_classified( ifg );

    % compute the core fiber from the fiber group (the tact profile is computed here)
    [fa, md, rd, ad, cl, core] = dtiComputeDiffusionPropertiesAlongFG( fg, dt,[],[],200);
    
    % How to make a trct profile from a NIFTI file (such as from a run model)
    % nifti_file = niftiRead('path/to/nifti/file.nii.gz')
    % val = dtiComputeDiffusionPropertiesAlongFG( fg, nifti_file,[],[],200);
    
    % 3. Select a center portion fo the tract and show the FA and MD values
    % normally we only use for analyses the middle most reliable portion of the fiber.
    nodesToPlot = 50:151;
    
    h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'off');
    tract_profile = plot(fa(nodesToPlot),'color', [0.2 0.2 0.9],'linewidth',4)
    set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
        'xticklabel',{'Tract begin','Tract end'},'xlim',[0 100],'ylim',[0.00 1.00],'Ytick',[0 .25 .5 .75],'Xtick',[0 100])
    Title_plot = title(fg.name);
    xlabel('Location on tract')
    ylh = ylabel('Fractional Anisotropy');
    
    saveas(tract_profile, fullfile('images/', Title_plot.String), 'png')
    saveas(tract_profile, fullfile('images/', Title_plot.String), 'eps')
    %clear()
    json.images(ifg).filename = strcat(Title_plot.String,'.png');
    json.images(ifg).name = strcat(Title_plot.String, ' tract analysis profile');
    json.images(ifg).desc = strcat(Title_plot.String, ' tract analysis profile');
    
end

savejson('', json, fullfile('images.json'));
end

