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

% make directories and set up variables
mkdir('images');
mkdir('profiles');
numfiles = 0;
imgnum = 0;
possible_error=0;
failed_tracts=[];


% load config.json
config = loadjson('config.json');

% load segmentation file and set number of nodes
load(fullfile(config.afq));
numnodes = config.numnodes;

% load tensor and noddi (if applicable) files
if isfield(config,'tensor')
    tensors = dir(fullfile(config.tensor,'*.nii.gz*'));
    tensors = [tensors(1) tensors(5) tensors(6) tensors(7)];
    tensors = tensors';
    end_index = 4;
else
    end_index = 0;
end

if isfield(config,'noddi')
    noddis = dir(fullfile(config.noddi,'*_NEW.nii.gz*'));
end

% Set data structures
if isfield(config,'tensor')
    for ii = 1:length(tensors)
        nii(ii).name = extractBefore(tensors(ii).name,strlength(tensors(ii).name)-6);
        nii(ii).data = niftiRead(fullfile(tensors(ii).folder,tensors(ii).name));
    end
end

if isfield(config,'noddi')
    for ii = 1:length(noddis)
        nii(4+ii).name = extractBetween(noddis(ii).name,'FIT_','_NEW');
        nii(4+ii).data = niftiRead(fullfile(noddis(ii).folder,noddis(ii).name));
    end
end

% Set up cell for csv
tract_profiles = cell(numnodes, length(nii));

for ifg = 1:length(fg_classified) 
try
    if config.fiberbased == 0
        display 'volume based statistics'
        fg = fg_classified( ifg );
        for jj = 1:length(nii)
            [tract, ~, ~, ~, ~, ~, ~, ~, ~, ~, myValsFgSTD] = dtiComputeDiffusionPropertiesAlongFG_sd( fg, nii(jj).data,[],[],numnodes);
            nii(jj).mean = tract;
            nii(jj).std = myValsFgSTD;
        end
    else
        display 'fiber based statistics'
        fgTract = fg_classified( ifg );
        fg = dtiXformFiberCoords(fgTract, inv(nii(end_index + 1).data.qto_xyz),'img'); % convert fibergroup to the proper space
        for jj = 1:length(nii)
            tract = Compute_FA_AlongFG(fg, nii(jj).data, [], [], numnodes);
            nii(jj).mean = nanmean(tract);
            nii(jj).std = nanstd(tract);
        end
    end
    
    for jj = 1:length(nii)
        tract_profiles(:,jj,1) = num2cell(nii(jj).mean);
        tract_profiles(:,jj,2) = num2cell(nii(jj).std);
    end
    
    for jj = 1:length(nii)
        T(:,jj) = table([tract_profiles(:,jj,1),tract_profiles(:,jj,2)]);
        T.Properties.VariableNames{jj} = char(nii(jj).name);
    end
    
    writetable(T, strcat('profiles/', strrep(fg.name, ' ', '_'), '_profiles.csv'));

    if config.ad == 1
        [imgnum,json] = analysisProfiles(nii(1).mean,fg,nii(1).name,'Axial Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,imgnum);
        numfiles = numfiles + 1;
    end
    
    if config.fa == 1
        [imgnum,json] = analysisProflies(nii(2).mean,fg,nii(2).name,'Fractional Anisotropy',[0.00, 1.00],[0 .25 .5 .75],imgnum);
        numfiles = numfiles + 1;
    end
    
    if config.md == 1
        [imgnum,json] = analysisProfiles(nii(3).mean,fg,nii(3).name,'Mean Diffusivity',[0.00, 2.00],[0 .5 1 1.5],imgnum);
        numfiles = numfiles + 1;
    end
    
    if config.rd == 1
        [imgnum,json] = analysisProfiles(nii(4).mean,fg,nii(4).name,'Radial Diffusivity',[0.00, 2.00],[0 .5 1 1.5],imgnum);
        numfiles = numfiles + 1;
    end
    
    if config.icvf == 1
        [imgnum,json] = analysisProfiles(nii(end_index + 1).mean,fg,nii(end_index + 1).name,'ICVF',[0 1.00],[0.25 .5 .75],imgnum);
        numfiles = numfiles + 1;
    end

    if config.isovf == 1
        [imgnum,json] = analysisProfiles(nii(end_index + 2).mean,fg,nii(end_index + 1).name,'ISOVF',[0 1.00],[0.25 .5 .75],imgnum);
        numfiles = numfiles + 1;
    end
    
    if config.od == 1
        [imgnum,json] = analysisProfiles(nii(end_index + 3).mean,fg,nii(end_index + 1).name,'OD',[0 1.00],[0.25 .5 .75],imgnum);
        numfiles = numfiles + 1;
    end
    
catch ME
    possible_error=1;
    failed_tracts = [failed_tracts, fg.name];
    
save('profiles/error_messages.mat','ME');
end
clf
end

fileID = fopen('numfiles.txt','w');
fprintf(fileID, '%d', numfiles);
fclose(fileID);

if possible_error==1
    results.quality_check = 'ERROR: The following tracts failed:';
    results.failed_tracts = failed_tracts;
else
    results.quality_check = 'All tracts analysis profiles were created successfully';
end
savejson('', results, 'product.json');

savejson('', json, fullfile('images.json'));
end

