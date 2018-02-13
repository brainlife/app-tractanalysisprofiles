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
noddi = {'ICVF','OD','ISOVF'};
noddiImage = {'ICVF','OD','ISOVF'};
nii_icvf = niftiRead(fullfile(config.icvf));
nii_od = niftiRead(fullfile(config.od));
nii_isovf = niftiRead(fullfile(config.isovf));
load(config.afq);
%load('output.mat');
numnodes = config.numnodes;

numfiles = 0;
if config.icvf
    numfiles = numfiles + length(fg_classified);
end
if config.od
    numfiles = numfiles + length(fg_classified);
end
if config.isovf
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
    [ICVF_tract, ICVF_SuperFiber, ~, ~] = Compute_FA_AlongFG(fg, nii_icvf, [], [], numnodes);
    mean_icvf = nanmean(ICVF_tract);
    
    [OD_tract, OD_SuperFiber, ~, ~] = Compute_FA_AlongFG(fg, nii_od, [], [], numnodes);
    mean_od = nanmean(OD_tract);
    
    [ISOVF_tract, ISOVF_SuperFiber, ~, ~] = Compute_FA_AlongFG(fg, nii_isovf, [], [], numnodes);
    mean_isovf = nanmean(ISOVF_tract);

    %[fa, md, rd, ad, cl, core] = dtiComputeDiffusionPropertiesAlongFG( fg, dt,[],[],100);
    tract_profiles = cell(numnodes, 3);
    

    tract_profiles(:,1) = num2cell(mean_icvf);
    tract_profiles(:,2) = num2cell(mean_od);
    tract_profiles(:,3) = num2cell(mean_isovf);
    
    T = cell2table(tract_profiles);
    T.Properties.VariableNames = {'ICVF','OD','ISOVF'};
    writetable(T, strcat('profiles/', strrep(fg.name, ' ', '_'), '_profiles.csv'));
    
    h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'off');
    
    if config.icvf
        tract_analysis_profiles(fg,mean_icvf,noddi{1},numnodes,'ICVF',[0 1.00],[0.25 .5 .75]);
        imgnum = imgnum+1;
        json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,noddi{1}),'.png');
        json.images(imgnum).name = strcat(fg.name);
        json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
        clf
    end
    if config.od
        tract_analysis_profiles(fg,mean_od,noddi{2},numnodes,'OD',[0 1.00],[0.25 .5 .75]);
        imgnum = imgnum+1;
        json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,noddi{2}),'.png');
        json.images(imgnum).name = strcat(fg.name);
        json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
        clf
    end
    if config.isovf
        tract_analysis_profiles(fg,mean_isovf,noddi{3},numnodes,'ISOVF',[0 1.00],[0.25 .5 .75]);
        imgnum = imgnum+1;
        json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,noddi{3}),'.png');
        json.images(imgnum).name = strcat(fg.name);
        json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
        clf
    end
    clf
    
catch ME
    possible_error=1;
    failed_tracts = [failed_tracts, fg.name];
    
save('profiles/error_messages.mat','ME')  
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

