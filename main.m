function [] = main()

if ~isdeployed
    disp('loading paths for IUHPC')
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/jsonlab'))

    disp('loading paths for Jetstream VM')
    addpath(genpath('/usr/local/vistasoft'))
    addpath(genpath('/usr/local/jsonlab'))
end

% make directories and set up variables
mkdir('images');
mkdir('profiles');
numfiles = 1;
possible_error=0;
failed_tracts=[];

% load config.json
config = loadjson('config.json');

if isempty(config.tensor) && isempty(config.noddi)
    display('No input initialized. Please specify input');
    exit
end

% load segmentation file and set number of nodes
load(fullfile(config.afq));
numnodes = config.numnodes;

% load tensor and noddi (if applicable) files
if isfield(config,'tensor')
    %tensors = dir(fullfile(config.tensor,'*.nii.gz*'));
    ad = dir(fullfile(config.tensor,'ad.nii.gz*'));
    fa = dir(fullfile(config.tensor,'fa.nii.gz*'));
    md = dir(fullfile(config.tensor,'md.nii.gz*'));
    rd = dir(fullfile(config.tensor,'rd.nii.gz*'));
    tensors = [ad(1) fa(1) md(1) rd(1)];
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
        nii(ii).name = char(extractBefore(tensors(ii).name,strlength(tensors(ii).name)-6));
        nii(ii).data = niftiRead(fullfile(tensors(ii).folder,tensors(ii).name));
        nii(ii).data_inv = 1./nii(ii).data.data;
        nii(ii).data_inv(~isfinite(nii(ii).data_inv))=0;
        if nii(ii).name == 'fa'
            nii(ii).units = 'unitless';
        else
            nii(ii).units = 'um^2/msec';
        end
        nii(end_index+ii).name = strcat(char(extractBefore(tensors(ii).name,strlength(tensors(ii).name)-6)),'_inverse');
        nii(end_index+ii).data = nii(ii).data;
        nii(end_index+ii).data.data = nii(ii).data_inv;
        nii(end_index+ii).data_inv = nii(ii).data_inv;
        if nii(end_index+ii).name == 'fa_inverse'
            nii(end_index+ii).units = 'unitless';
        else
            nii(end_index+ii).units = 'msec/um^2';
        end
    end
    end_index = length(nii);
end

if isfield(config,'noddi')
    for ii = 1:length(noddis)
        nii(end_index+ii).name = char(extractBetween(noddis(ii).name,'FIT_','_NEW'));
        nii(end_index+ii).data = niftiRead(fullfile(noddis(ii).folder,noddis(ii).name));
        nii(end_index+ii).data_inv = 1./nii(ii).data.data;
        nii(end_index+ii).data_inv(~isfinite(nii(end_index+ii).data_inv))=0;
        nii(end_index+ii).units = 'unitless';
    end
    end_index = length(nii);
    for ii = 1:length(noddis)
        nii(end_index+ii).name = strcat(char(extractBetween(noddis(ii).name,'FIT_','_NEW')),'_inverse');
        nii(end_index+ii).data = nii(end_index-3+ii).data;
        nii(end_index+ii).data_inv = nii(end_index-3+ii).data_inv;
        nii(end_index+ii).data.data = nii(end_index+ii).data_inv;
        nii(end_index+ii).units = 'unitless';
    end
    end_index = length(nii);
end

% Set up cell for csv
tract_profiles = cell(numnodes, length(nii));

for ifg = 1:length(fg_classified) 
    try
        if config.fiberbased == 0
            display 'volume based statistics'
            fg = fg_classified( ifg );
            for jj = 1:length(nii)
                display(sprintf('computing %s',nii(jj).name));
                [tract, ~, ~, ~, ~, ~, ~, ~, ~, ~, myValsFgSTD] = dtiComputeDiffusionPropertiesAlongFG_sd( fg, nii(jj).data,[],[],numnodes);
                nii(jj).mean = tract;
                nii(jj).std = myValsFgSTD;
            end
        else
            display 'fiber based statistics'
            fgTract = fg_classified( ifg );
            fg = dtiXformFiberCoords(fgTract, inv(nii(2).data.qto_xyz),'img'); % convert fibergroup to the proper space
            for jj = 1:length(nii)
                display(sprintf('computing %s',nii(jj).name));
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
            T.Properties.VariableUnits{jj} = nii(jj).units;
        end
        
        writetable(T, strcat('profiles/', strrep(fg.name, ' ', '_'), '_profiles.csv'));

        if isfield(config,'tensor')
            % AD
            analysisProfiles(nii(1).mean,fg,nii(1).name,'Axial Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(1).units);
            json.images(numfiles).filename = strcat('images/',fg.name,'_', 'Axial Diffusivity','.png');
            json.images(numfiles).name = strcat(fg.name);
            json.images(numfiles).desc = strcat(fg.name, ' tract analysis profile');
            numfiles = numfiles + 1;
            % FA
            analysisProfiles(nii(2).mean,fg,nii(2).name,'Fractional Anisotropy',[0.00, 1.00],[0 .25 .5 .75],numnodes,nii(2).units);
            json.images(numfiles).filename = strcat('images/',fg.name,'_', 'Fractional Anisotropy','.png');
            json.images(numfiles).name = strcat(fg.name);
            json.images(numfiles).desc = strcat(fg.name, ' tract analysis profile');
            numfiles = numfiles + 1;
            % MD
            analysisProfiles(nii(3).mean,fg,nii(3).name,'Mean Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(3).units);
            json.images(numfiles).filename = strcat('images/',fg.name,'_', 'Mean Diffusivity','.png');
            json.images(numfiles).name = strcat(fg.name);
            json.images(numfiles).desc = strcat(fg.name, ' tract analysis profile');
            numfiles = numfiles + 1;
            % RD
            analysisProfiles(nii(4).mean,fg,nii(4).name,'Radial Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(4).units);
            json.images(numfiles).filename = strcat('images/',fg.name,'_', 'Radial Diffusivity','.png');
            json.images(numfiles).name = strcat(fg.name);
            json.images(numfiles).desc = strcat(fg.name, ' tract analysis profile');
            numfiles = numfiles + 1;
        end
        
        if isfield(config,'noddi')
            % ICVF
            analysisProfiles(nii(end_index-6+1).mean,fg,nii(end_index-6+1).name,'ICVF',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+1).units);
            json.images(numfiles).filename = strcat('images/',fg.name,'_', 'ICVF','.png');
            json.images(numfiles).name = strcat(fg.name);
            json.images(numfiles).desc = strcat(fg.name, ' tract analysis profile');
            numfiles = numfiles + 1;
            % ISOVF
            analysisProfiles(nii(end_index-6+2).mean,fg,nii(end_index-6+2).name,'ISOVF',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+2).units);
            json.images(numfiles).filename = strcat('images/',fg.name,'_', 'ISOVF','.png');
            json.images(numfiles).name = strcat(fg.name);
            json.images(numfiles).desc = strcat(fg.name, ' tract analysis profile');
            numfiles = numfiles + 1;
            % OD
            analysisProfiles(nii(end_index-6+3).mean,fg,nii(end_index-6+3).name,'OD',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+3).units);
            json.images(numfiles).filename = strcat('images/',fg.name,'_', 'OD','.png');
            json.images(numfiles).name = strcat(fg.name);
            json.images(numfiles).desc = strcat(fg.name, ' tract analysis profile');
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
fprintf(fileID, '%d', numfiles-1); %matlab uses 1 based indexing
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

