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

imgnum = 0;

possible_error=0;
failed_tracts=[];

% load my own config.json
config = loadjson('config.json');
tensor = {'FA','MD','RD','AD'};
noddi = {'ICVF','OD','ISOVF'};
noddiImage = {'ICVF','OD','ISOVF'};
mkdir('images');
mkdir('profiles');
load(config.afq);

numnodes = config.numnodes;

numfiles = 0;

if isfield(config,'dt6') ~= 1
    display('No Tensor.')
else
    if config.fa
        dt = dtiLoadDt6(fullfile(config.dt6,'/dti/dt6.mat'));
        numfiles = numfiles + length(fg_classified);
    end

    if config.md
        dt = dtiLoadDt6(fullfile(config.dt6,'/dti/dt6.mat'));
        numfiles = numfiles + length(fg_classified);
    end

    if config.rd
        dt = dtiLoadDt6(fullfile(config.dt6,'/dti/dt6.mat'));
        numfiles = numfiles + length(fg_classified);
    end

    if config.ad
        dt = dtiLoadDt6(fullfile(config.dt6,'/dti/dt6.mat'));
        numfiles = numfiles + length(fg_classified);
    end
end
    
if isfield(config,'noddi') ~= 1
    display('No NODDI')
else
    if config.icvf
        noddiImage{1} = niftiRead(fullfile(config.noddi,'FIT_ICVF_NEW.nii.gz'));
        numfiles = numfiles + length(fg_classified);
    end

    if config.od
        noddiImage{2} = niftiRead(fullfile(config.noddi,'FIT_OD_NEW.nii.gz'));
        numfiles = numfiles + length(fg_classified);
    end

    if config.isovf
        noddiImage{3} = niftiRead(fullfile(config.noddi,'FIT_ISOVF_NEW.nii.gz'));
        numfiles = numfiles + length(fg_classified);
    end
end

fileID = fopen('numfiles.txt','w');
fprintf(fileID, '%d', numfiles);
fclose(fileID);

for ifg = 1:length(fg_classified)
try
    fg = fg_classified( ifg );

	if isfield(config,'dt6') == 1
		% compute the core fiber from the fiber group (the tact profile is computed here)
		[fa, md, rd, ad, cl, SuperFiber, fgClipped, cp, cs, fgResampled] = dtiComputeDiffusionPropertiesAlongFG( fg, dt,[],[],numnodes);
    		
		tract_profiles = cell(numnodes, 4);
    

		tract_profiles(:,1) = num2cell(fa);
		tract_profiles(:,2) = num2cell(md);
		tract_profiles(:,3) = num2cell(rd);
		tract_profiles(:,4) = num2cell(ad);
    
		T = cell2table(tract_profiles);
		T.Properties.VariableNames = tensor;
		writetable(T, strcat('profiles/', strrep(fg.name, ' ', '_'), '_tensor_profiles.csv'));
		if config.fa
			tract_analysis_profiles(fg,fa,tensor{1},numnodes,'Fractional Anisotropy',[0 1.00],[0 .25 .5 .75]);
            		imgnum = imgnum+1;
            		json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,tensor{1}),'.png');
            		json.images(imgnum).name = strcat(fg.name);
            		json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
            		clf
		end
		if config.md
			tract_analysis_profiles(fg,md,tensor{2},numnodes,'Mean Diffusivity',[0 1.00],[0 .25 .5 .75]);
            		imgnum = imgnum+1;
            		json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,tensor{2}),'.png');
            		json.images(imgnum).name = strcat(fg.name);
            		json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
            		clf		
		end
		if config.rd
			tract_analysis_profiles(fg,rd,tensor{3},numnodes,'Radial Diffusivity',[0 1.00],[0 .25 .5 .75]);
            		imgnum = imgnum+1;
            		json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,tensor{3}),'.png');
            		json.images(imgnum).name = strcat(fg.name);
            		json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
            		clf
		end
		if config.ad
			tract_analysis_profiles(fg,ad,tensor{4},numnodes,'Axial Diffusivity',[0 2.00],[0 .5 1 1.5]);
            		imgnum = imgnum+1;
            		json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,tensor{4}),'.png');
            		json.images(imgnum).name = strcat(fg.name);
            		json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
            		clf
        	end
	end

	if isfield(config,'noddi') == 1
        	% compute the core fiber from the fiber group (the tact profile is computed here)
        	[fa1, SuperFiber1, fgClipped1, cp1, cs1, fgResampled1] = dtiComputeDiffusionPropertiesAlongFG( fg, noddiImage{1},[],[],numnodes);
        	[fa2, SuperFiber2, fgClipped2, cp2, cs2, fgResampled2] = dtiComputeDiffusionPropertiesAlongFG( fg, noddiImage{2},[],[],numnodes);
        	[fa3, SuperFiber3, fgClipped3, cp3, cs3, fgResampled3] = dtiComputeDiffusionPropertiesAlongFG( fg, noddiImage{3},[],[],numnodes);

        	tract_profiles = cell(numnodes, 3);
    

        	tract_profiles(:,1) = num2cell(fa1);
        	tract_profiles(:,2) = num2cell(fa2);
        	tract_profiles(:,3) = num2cell(fa3);
    
        	T = cell2table(tract_profiles);
        	T.Properties.VariableNames = noddi;
        	writetable(T, strcat('profiles/', strrep(fg.name, ' ', '_'), '_noddi_profiles.csv'));
        	if config.icvf
            		tract_analysis_profiles(fg,fa1,noddi{1},numnodes,'ICVF',[0 1.00],[0.25 .5 .75]);
            		imgnum = imgnum+1;
            		json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,noddi{1}),'.png');
            		json.images(imgnum).name = strcat(fg.name);
            		json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
            		clf
		end
        	if config.icvf
            		tract_analysis_profiles(fg,fa2,noddi{2},numnodes,'OD',[0 1.00],[0.25 .5 .75]);
            		imgnum = imgnum+1;
            		json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,noddi{2}),'.png');
            		json.images(imgnum).name = strcat(fg.name);
            		json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
            		clf
        	end
        	if config.icvf
            		tract_analysis_profiles(fg,fa3,noddi{3},numnodes,'ISOVF',[0 1.00],[0.25 .5 .75]);
            		imgnum = imgnum+1;
            		json.images(imgnum).filename = strcat('images/',sprintf('%s_%s',fg.name,noddi{3}),'.png');
            		json.images(imgnum).name = strcat(fg.name);
            		json.images(imgnum).desc = strcat(fg.name, ' tract analysis profile');
            		clf
        	end
        	clf
    	end

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
