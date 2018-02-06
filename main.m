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
%dt = dtiLoadDt6(fullfile(config.dt6,'/dti/dt6.mat'));
nii_fa = niftiRead(fullfile(config.fa));
nii_md = niftiRead(fullfile(config.md));
nii_rd = niftiRead(fullfile(config.rd));
nii_ad = niftiRead(fullfile(config.ad));
load(config.afq);
%load('output.mat');
numnodes = config.numnodes;

numfiles = 0;
if config.fa
    numfiles = numfiles + length(fg_classified);
end
if config.md
    numfiles = numfiles + length(fg_classified);
end
if config.rd
    numfiles = numfiles + length(fg_classified);
end
if config.ad
    numfiles = numfiles + length(fg_classified);
end

fileID = fopen('numfiles.txt','w');
fprintf(fileID, '%d', numfiles);
fclose(fileID);

mkdir('images');
mkdir('profiles');
imgnum = 0;

possible_error=0;
failed_tracts=[];
for ifg = 1:length(fg_classified)
try
    fgTract = fg_classified( ifg );
    fg = dtiXformFiberCoords(fgTract, inv(nii_fa.qto_xyz),'img'); % convert fibergroup to the proper space
    
    % compute the core fiber from the fiber group (the tact profile is computed here)
    [FA_tract, FA_SuperFiber, ~, ~] = Compute_FA_AlongFG(fg, nii_fa, [], [], numnodes);
    mean_fa = nanmean(FA_tract);
    
    [MD_tract, MD_SuperFiber, ~, ~] = Compute_FA_AlongFG(fg, nii_md, [], [], numnodes);
    mean_md = nanmean(MD_tract);
    
    [RD_tract, RD_SuperFiber, ~, ~] = Compute_FA_AlongFG(fg, nii_rd, [], [], numnodes);
    mean_rd = nanmean(RD_tract);
    
    [AD_tract, AD_SuperFiber, ~, ~] = Compute_FA_AlongFG(fg, nii_ad, [], [], numnodes);
    mean_ad = nanmean(AD_tract);

    %[fa, md, rd, ad, cl, core] = dtiComputeDiffusionPropertiesAlongFG( fg, dt,[],[],100);
    tract_profiles = cell(numnodes, 4);
    

    tract_profiles(:,1) = num2cell(mean_fa);
    tract_profiles(:,2) = num2cell(mean_md);
    tract_profiles(:,3) = num2cell(mean_rd);
    tract_profiles(:,4) = num2cell(mean_ad);
    
    T = cell2table(tract_profiles);
    T.Properties.VariableNames = {'FA', 'MD', 'RD', 'AD'};
    writetable(T, strcat('profiles/', strrep(fg.name, ' ', '_'), '_profiles.csv'));
    % How to make a trct profile from a NIFTI file (such as from a run model)
    % nifti_file = niftiRead('path/to/nifti/file.nii.gz')
    % val = dtiComputeDiffusionPropertiesAlongFG( fg, nifti_file,[],[],200);
    
    % 3. Select a center portion fo the tract and show the FA and MD values
    % normally we only use for analyses the middle most reliable portion of the fiber.
    %nodesToPlot = 1:200;
    %nodesToPlot = 25:76;
    
    h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'off');
    

    if config.fa
        tract_profile = plot(mean_fa,'color', [0.2 0.2 0.9],'linewidth',4)
        ylh = ylabel('Fractional Anisotropy');
        ylim = [0.00, 1.00];
        ytick = [0 .25 .5 .75];
        set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
        'xticklabel',{'Tract begin','Tract end'},'xlim',[0 numnodes],'ylim',ylim,'Ytick',ytick,'Xtick',[0 numnodes])
        Title_plot = title(fg.name);
        xlabel('Location on tract')
        saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_fa')), 'png')
        saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_fa')), 'eps')
        %clear()
        imgnum = imgnum+1;
        json.images(imgnum).filename = strcat('images/',Title_plot.String,'_fa','.png');
        json.images(imgnum).name = strcat(Title_plot.String);
        json.images(imgnum).desc = strcat(Title_plot.String, ' tract analysis profile');
        clf
    end
    if config.md
        tract_profile = plot(mean_md,'color', [0.2 0.2 0.9],'linewidth',4)
        ylh = ylabel('Mean Diffusivity');
        ylim = [0.00, 2.00];
        ytick = [0 .5 1 1.5];
        set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
        'xticklabel',{'Tract begin','Tract end'},'xlim',[0 numnodes],'ylim',ylim,'Ytick',ytick,'Xtick',[0 numnodes])
        Title_plot = title(fg.name);
        xlabel('Location on tract')
        saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_md')), 'png')
        saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_md')), 'eps')
        %clear()
        imgnum = imgnum + 1;
        json.images(imgnum).filename = strcat('images/',Title_plot.String,'_md', '.png');
        json.images(imgnum).name = strcat(Title_plot.String);
        json.images(imgnum).desc = strcat(Title_plot.String, ' tract analysis profile');
        clf
    end
    if config.rd
        tract_profile = plot(mean_rd,'color', [0.2 0.2 0.9],'linewidth',4)
        ylh = ylabel('Radial Diffusivity');
        ylim = [0.00, 2.00];
        ytick = [0 .5 1 1.5];
        set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
        'xticklabel',{'Tract begin','Tract end'},'xlim',[0 numnodes],'ylim',ylim,'Ytick',ytick,'Xtick',[0 numnodes])
        Title_plot = title(fg.name);
        xlabel('Location on tract')
        saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_rd')), 'png')
        saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_rd')), 'eps')
        %clear()
        imgnum = imgnum + 1;
        json.images(imgnum).filename = strcat('images/',Title_plot.String,'_rd', '.png');
        json.images(imgnum).name = strcat(Title_plot.String);
        json.images(imgnum).desc = strcat(Title_plot.String, ' tract analysis profile');
        clf
    end
    if config.ad
        tract_profile = plot(mean_ad,'color', [0.2 0.2 0.9],'linewidth',4)
        ylh = ylabel('Axial Diffusivity');
        ylim = [0.00, 2.00];
        ytick = [0 .5 1 1.5];
        set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
        'xticklabel',{'Tract begin','Tract end'},'xlim',[0 numnodes],'ylim',ylim,'Ytick',ytick,'Xtick',[0 numnodes])
        Title_plot = title(fg.name);
        xlabel('Location on tract')
        saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_ad')), 'png')
        saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_ad')), 'eps')
        %clear()
        imgnum = imgnum + 1;
        json.images(imgnum).filename = strcat('images/',Title_plot.String,'_ad', '.png');
        json.images(imgnum).name = strcat(Title_plot.String);
        json.images(imgnum).desc = strcat(Title_plot.String, ' tract analysis profile');
        clf
    end
catch ME
    possible_error=1;
    failed_tracts = [failed_tracts, fg.name];
    
save('profiles/error_messages.mat','ME')

%     set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
%         'xticklabel',{'Tract begin','Tract end'},'xlim',[0 50],'ylim',ylim,'Ytick',ytick,'Xtick',[0 50])
%     Title_plot = title(fg.name);
%     xlabel('Location on tract')
%     
%     
%     saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_', prop)), 'png')
%     saveas(tract_profile, fullfile('images/', strcat(Title_plot.String, '_', prop)), 'eps')
%     %clear()
%     json.images(ifg).filename = strcat('images/',Title_plot.String,'_', prop, '.png');
%     json.images(ifg).name = strcat(Title_plot.String);
%     json.images(ifg).desc = strcat(Title_plot.String, ' tract analysis profile');
%     
end
end
clf

if possible_error==1
    results.quality_check = 'ERROR: The following tracts failed:';
    results.failed_tracts = failed_tracts;
else
    results.quality_check = 'All tracts analysis profiles were created successfully';
end
savejson('', results, 'product.json');

savejson('', json, fullfile('images.json'));
end

