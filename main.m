function [] = main()

if ~isdeployed
    disp('loading paths for IUHPC')
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/soft/mason/SPM/spm8'))
    addpath(genpath('/N/u/brlife/git/jsonlab'))

    disp('loading paths for Jetstream VM')
    addpath(genpath('/usr/local/vistasoft'))
    addpath(genpath('/usr/local/spm8'))
    addpath(genpath('/usr/local/jsonlab'))
end

% make directories and set up variables
mkdir('images');
mkdir('profiles');
numfiles = 1;
possible_error=0;
possible_error_lows=0;
failed_tracts=[];
failed_tracts_lows=[];

% load config.json
config = loadjson('config.json');

if ~isfield(config,'ad') && ~isfield(config,'icvf')
    display('No input initialized. Please specify input');
    exit
end

% load segmentation file and set number of nodes
load(fullfile(config.afq));
numnodes = config.numnodes;

if ~exist('fg_classified','var')
    fg_classified = {tracts};
elseif ~iscell(fg_classified)
    fg_classified = {fg_classified};
else
    fg_classified = fg_classified;
end

% load tensor and noddi (if applicable) files
if isfield(config,'ad')
    %tensors = dir(fullfile(config.tensor,'*.nii.gz*'));
    ad = dir(config.ad);
    fa = dir(config.fa);
    md = dir(config.md);
    rd = dir(config.rd);
    tensors = [ad fa md rd];
    end_index = 4;
else
    ad = [];
    end_index = 0;
end

if isfield(config,'icvf')
    icvf = dir(config.icvf);
    isovf = dir(config.isovf);
    od = dir(config.od);
    noddis = [icvf isovf od];
end

% Set data structures
if isfield(config,'ad')
    for ii = 1:length(tensors)
        nii(ii).name = char(extractBefore(tensors(ii).name,strlength(tensors(ii).name)-6));
        nii(ii).data = niftiRead(fullfile(tensors(ii).folder,tensors(ii).name));
        nii(ii).non_zero_index = find(nii(ii).data.data(:,:,:) ~= 0);
        if max(nii(ii).data.data(nii(ii).non_zero_index)) < 0.01 && ~strcmp(nii(ii).name,'fa')
            nii(ii).data.data = nii(ii).data.data * 1000;
        end
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

if isfield(config,'icvf')
    for ii = 1:length(noddis)
        nii(end_index+ii).name = char(extractBefore(noddis(ii).name,strlength(noddis(ii).name)-6));
        nii(end_index+ii).data = niftiRead(fullfile(noddis(ii).folder,noddis(ii).name));
        nii(end_index+ii).data_inv = 1./nii(ii).data.data;
        nii(end_index+ii).data_inv(~isfinite(nii(end_index+ii).data_inv))=0;
        nii(end_index+ii).units = 'unitless';
    end
    end_index = length(nii);
    for ii = 1:length(noddis)
        nii(end_index+ii).name = strcat(char(extractBefore(noddis(ii).name,strlength(noddis(ii).name)-6)),'_inverse');
        nii(end_index+ii).data = nii(end_index-3+ii).data;
        nii(end_index+ii).data_inv = nii(end_index-3+ii).data_inv;
        nii(end_index+ii).data.data = nii(end_index+ii).data_inv;
        nii(end_index+ii).units = 'unitless';
    end
    end_index = length(nii);
end

% Set up cell for csv
tract_profiles = cell(numnodes, length(nii));

for ifg = 1:length(fg_classified{1})
    try
        if config.fiberbased == 0
            display 'volume based statistics'
            [SuperFiber, fgResampled] = dtiComputeSuperFiberRepresentation(fg_classified{1}(ifg), [], numnodes);
            for jj = 1:length(nii)
                if length(fgResampled.fibers) < 6
                    display('too few streamlines. outputting profile of NaNs')
                    nii(jj).mean = NaN(numnodes,1);
                    nii(jj).std = NaN(numnodes,1);
                else
                    display(sprintf('computing %s',nii(jj).name));
                    [tract, ~, ~, ~, ~, ~, ~, ~, ~, ~, myValsFgSTD] = dtiComputeDiffusionPropertiesAlongFG_sd( fgResampled, nii(jj).data,[],[],numnodes,[],SuperFiber);
                    nii(jj).mean = tract;
                    nii(jj).std = myValsFgSTD;
                end
                clear tract myValsFgSTD
            end
        else
            display 'fiber based statistics'
            [SuperFiber, fgResampled] = dtiComputeSuperFiberRepresentation(dtiXformFiberCoords(fg_classified{1}(ifg), inv(nii(2).data.qto_xyz),'img'),[], numberOfNodes); % convert fibergroup to the proper space
            for jj = 1:length(nii)
                if length(fgResampled.fibers) < 6
                    display('too few streamlines. outputting profile of NaNs')
                    nii(jj).mean = NaN(numnodes,1);
                    nii(jj).std = NaN(numnodes,1);
                else
                    display(sprintf('computing %s',nii(jj).name));
                    tract = Compute_FA_AlongFG(fgResampled, nii(jj).data, [], [], numnodes,SuperFiber);
                    nii(jj).mean = nanmean(tract);
                    nii(jj).std = nanstd(tract);
                end
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
        
        fg_filename = strrep(fgResampled.name, ' ', '_');
        writetable(T, strcat('profiles/', fg_filename, '_profiles.csv'));

        if isfield(config,'ad')
            % AD
            analysisProfiles(nii(1).mean,fgResampled,nii(1).name,'Axial Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(1).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_ad.png');
            json.images(numfiles).name = fgResampled.name;
            json.images(numfiles).desc = strcat('Axial Diffusivity');
            numfiles = numfiles + 1;
            % FA
            analysisProfiles(nii(2).mean,fgResampled,nii(2).name,'Fractional Anisotropy',[0.00, 1.00],[0 .25 .5 .75],numnodes,nii(2).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_fa.png');
            json.images(numfiles).name = fgResampled.name;
            json.images(numfiles).desc = strcat('Fractional Anistropy');
            numfiles = numfiles + 1;
            % MD
            analysisProfiles(nii(3).mean,fgResampled,nii(3).name,'Mean Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(3).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_md.png');
            json.images(numfiles).name = fgResampled.name;
            json.images(numfiles).desc = strcat('Mean Diffusivity');
            numfiles = numfiles + 1;
            % RD
            analysisProfiles(nii(4).mean,fgResampled,nii(4).name,'Radial Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(4).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_rd.png');
            json.images(numfiles).name = fgResampled.name;
            json.images(numfiles).desc = strcat('Radial Diffusivity');
            numfiles = numfiles + 1;
        end
        
        if isfield(config,'icvf')
            % ICVF
            analysisProfiles(nii(end_index-6+1).mean,fgResampled,nii(end_index-6+1).name,'ICVF',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+1).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_icvf.png');
            json.images(numfiles).name = fgResampled.name;
            json.images(numfiles).desc = strcat('ICVF');
            numfiles = numfiles + 1;
            % ISOVF
            analysisProfiles(nii(end_index-6+2).mean,fgResampled,nii(end_index-6+2).name,'ISOVF',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+2).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_isovf.png');
            json.images(numfiles).name = fgResampled.name;
            json.images(numfiles).desc = strcat('ISOVF');
            numfiles = numfiles + 1;
            % OD
            analysisProfiles(nii(end_index-6+3).mean,fgResampled,nii(end_index-6+3).name,'OD',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+3).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_od.png');
            json.images(numfiles).name = fgResampled.name;
            json.images(numfiles).desc = strcat('OD');
            numfiles = numfiles + 1;
        end
        
    catch ME
        possible_error=1;
        failed_tracts = [failed_tracts, fgResampled.name];
        
    save('profiles/error_messages.mat','ME');
    end
    
    if length(fgResampled.fibers) < 6
        possible_error_lows=1;
        failed_tracts_lows = [failed_tracts, fgResampled.name];
        save('profiles/error_messages_lows.mat','failed_tracts_lows')
    end
    
    clf
    clear fgResampled SuperFiber;
end

fileID = fopen('numfiles.txt','w');
fprintf(fileID, '%d', numfiles-1); %matlab uses 1 based indexing
fclose(fileID);

if possible_error==1
    results.quality_check = 'ERROR: The following tracts failed:';
    results.failed_tracts = failed_tracts;
elseif possible_error_lows==1
    results.quality_check = 'ERROR: The following tracts have too few streamlines:';
    results.failed_tracts = failed_tracts_lows;
else
    results.quality_check = 'All tracts analysis profiles were created successfully';
end


savejson('', results, 'product.json');
savejson('', json, fullfile('images.json'));

end

