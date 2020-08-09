function [] = main()

if ~isdeployed
    disp('loading paths for IUHPC')
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/soft/mason/SPM/spm8'))
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/soft/rhel7/mrtrix/3.0/mrtrix3/matlab'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
end

% make directories and set up variables
mkdir('images');
mkdir('profiles');
tensor = {'ad','fa','md','rd'};
dki = {'ga','mk','ak','rk'};
noddi = {'ndi','isovf','odi'};

% load config.json
config = loadjson('config.json');

if ~isfield(config,'ad') && ~isfield(config,'ndi')
    display('Please specify either tensor or noddi input (or both). You are trying to run this app with neither of them.');
    exit
end

% load segmentation file and set number of nodes; take in both
% classification structure and tck tractogram, will generate fg_classified
% structure
load(fullfile(config.afq));
numnodes = config.numnodes;
wbFG = wma_loadTck(config.tck);
fg_classified = bsc_makeFGsFromClassification_v5(classification, wbFG);

%%%% load tensor and noddi (if applicable) files
measures = {};
end_index = [];
scale_index = [];
value_units = [];
inverse_units = [];

% dti
if isfield(config,tensor(1))
	for ii = 1:length(tensor)
		measures{ii} = dir(config.(tensor{ii}));
    end
    end_index = [length(tensor)];
    scale_index = ["true","false","true","true"];
    value_units = ["um^2/msec","unitless","um^2/msec","um^2/msec"];
    inverse_units = ["msec/um^2","unitless","msec/um^2","um^2/msec"];

    % dki
    if isfield(config,dki(1))
        for kk = 1:length(dki)
            measures{ii+kk} = dir(config.(dki{kk}));
        end
        end_index = [end_index end_index+length(dki)];
        scale_index = [scale_index ["false","false","false","false"]];
        value_units = [value_units ["unitless","um^2/msec","um^2/msec","um^2/msec"]];
    	inverse_units = [inverse_units ["unitless","msec/um^2","msec/um^2","um^2/msec"]];
    end
else
    end_index = 0;
end

% noddi
if isfield(config,noddi(1))
	for nn = 1:length(noddi)
		measures{end_index(end)+nn} = dir(config.(noddi{nn}));
    end
end_index = [end_index end_index(end)+length(noddi)];
scale_index = [scale_index ["false","false","false"]];
value_units = [value_units ["unitless","unitless","unitless"]];
inverse_units = [inverse_units ["unitless","unitless","unitless"]];
end

%%%% load nifti data into singular structure
nii = build_nifti_data(measures,scale_index,value_units,inverse_units);

%%%% set up array for product.json
for ii = 1:length(classification.names)
    tractname = strrep(strrep(classification.names{ii},'.','_'),' ','');
    tractprofiles.(tractname) = struct();
    for jj = 1:length(nii)
        if ~strcmp(nii(jj).name(end),'e')
            measurename = nii(jj).name;
            tractprofiles.(tractname).(measurename).profile = [];
            tractprofiles.(tractname).(measurename).mean = [];
            tractprofiles.(tractname).(measurename).sd = [];
        end
    end
end

%%%% generate analysis profiles
tractAnalysisProfiles(classification,nii,config,numnodes,fg_classified)

end

