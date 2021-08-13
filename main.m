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
myelin = {'myelin'};
qmri = {'T1','R1','M0','PD','MTV','VIP','SIR','WF'};

% load config.json
config = loadjson('config.json');

if ~isfield(config,'ad') && ~isfield(config,'ndi') && ~isfield(config,'myelin')
    display('Please specify either tensor, noddi or myelin input (or both). You are trying to run this app with neither of them.');
    exit
end

% load segmentation file and set number of nodes; take in both
% classification structure and tck tractogram
load(fullfile(config.afq));
numnodes = config.numnodes;
wbFG = wma_loadTck(config.tck);

% identify if classification structure contains empty streamlines
indices = unique(classification.index(classification.index>0));
empty_indices = [];
j=1;
for i = 1:length(classification.names)
    if ~ismissing(i,indices)
        empty_indices(j) = i;
        j=j+1;
    end
end

% clean up before creation of fg_classified structure. will fail if empty
% classifications are included
cleaned_classification = classification;
if ~isempty(empty_indices)
    cleaned_classification.names = {classification.names{indices}};
    new_indices = [1:1:length(cleaned_classification.names)];
    for i = 1:length(indices)
        tract_indices = find(cleaned_classification.index==indices(i));
        for j = 1:length(tract_indices)
            cleaned_classification.index(tract_indices(j)) = new_indices(i);
        end
    end
end

fg_classified = bsc_makeFGsFromClassification_v5(cleaned_classification, wbFG);

%%%% load tensor and noddi (if applicable) files
measures = {};
end_index = [];
scale_index = [];
value_units = [];
inverse_units = [];

% dti
if exist(config.(tensor{1}))
	for ii = 1:length(tensor)
		measures{ii} = dir(config.(tensor{ii}));
    end
    end_index = [length(tensor)];
    scale_index = ["true","false","true","true"];
    value_units = ["um^2/msec","unitless","um^2/msec","um^2/msec"];
    inverse_units = ["msec/um^2","unitless","msec/um^2","um^2/msec"];

    % dki
    if exist(config.(dki{1}))
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

% myelin
if isfield(config,myelin(1))
    measures{end_index(end)+1} = dir(config.myelin);
    end_index = [end_index end_index(end)+1];
    scale_index = [scale_index ["false"]];
    value_units = [value_units ["unitless"]];
    inverse_units = [inverse_units ["unitless"]];
end

% qmri
if isfield(config,qmri{1})
    n=0
    for nn = 1:length(qmri)
        if exist(qmri{nn},'file')
            n=n+1
            measures{end_index(end)+n} = dir(config.(qmri{nn}));

            if strcmp(qmri{nn},'T1')
                scale_index = [scale_index ["false"]];
                value_units = [values_units ["seconds"]];
                inverse_units = [inverse_units ["1/seconds"]];
            elseif strcmp(qmri{nn},'R1')
               scale_index = [scale_index ["false"]];
                value_units = [values_units ["1/seconds"]];
                inverse_units = [inverse_units ["seconds"]];
            elseif strcmp(qmri{nn},'M0')
               scale_index = [scale_index ["false"]];
                value_units = [values_units ["unitless"]];
                inverse_units = [inverse_units ["unitless"]];
            elseif strcmp(qmri{nn},'MTV')
               scale_index = [scale_index ["false"]];
                value_units = [values_units ["unitless"]];
                inverse_units = [inverse_units ["unitless"]];
            elseif strcmp(qmri{nn},'VIP')
               scale_index = [scale_index ["false"]];
                value_units = [values_units ["unitless"]];
                inverse_units = [inverse_units ["unitless"]];
            elseif strcmp(qmri{nn},'SIR')
               scale_index = [scale_index ["false"]];
                value_units = [values_units ["unitless"]];
                inverse_units = [inverse_units ["unitless"]];
            elseif strcmp(qmri{nn},'WF')
               scale_index = [scale_index ["false"]];
                value_units = [values_units ["unitless"]];
                inverse_units = [inverse_units ["unitless"]];
            end
        end
    end
    end_index = [end_index end_index(end)+n];
end

display(measures{length(measures)})

%%%% load nifti data into singular structure
nii = build_nifti_data(measures,scale_index,value_units,inverse_units);

%%%% generate analysis profiles
tractAnalysisProfiles(classification,nii,config,numnodes,fg_classified,empty_indices)

end

