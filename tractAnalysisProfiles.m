function [] = tractAnalysisProfiles(classification,nii,config,numnodes,fg_classified,empty_indices)

numfiles = 1;
possible_error=0;
possible_error_lows=0;
failed_tracts=[];
failed_tracts_lows=[];
end_index=length(nii);
tps = {'mean','sd'};

for ii = 1:length(classification.names)
    tractname = strrep(strrep(classification.names{ii},'.','_'),' ','');
    tractprofiles.(tractname) = struct();
    for jj = 1:length(nii)
        if ismember(ii,empty_indices)
            % capture cases where classification structure has no streamlines
            for tp = 1:length(tps)
                T(:,(2*jj-2)+tp) = table(NaN(numnodes,1));
                T.Properties.VariableNames{(2*jj-2)+tp} = sprintf('%s_%s',char(nii(jj).name),tps{tp});
                T.Properties.VariableUnits{(2*jj-2)+tp} = sprintf('%s',nii(jj).units);
            end
            T.x_coords = NaN(numnodes,1);
            T.Properties.VariableUnits{length(T.Properties.VariableNames)} = 'mm';
            T.y_coords = NaN(numnodes,1);
            T.Properties.VariableUnits{length(T.Properties.VariableNames)} = 'mm';
            T.z_coords = NaN(numnodes,1);
            T.Properties.VariableUnits{length(T.Properties.VariableNames)} = 'mm';
            writetable(T, strcat('profiles/', tractname, '_profiles.csv'));
            
            possible_error_lows=1;
            failed_tracts_lows = strcat(failed_tracts,tractname," ");
            save('profiles/error_messages_lows.mat','failed_tracts_lows')
        end
        % generate tractprofiles structure, outputting nans for tracts with
        % zero streamlines
        if ~strcmp(nii(jj).name(end),'e')
            if ismember(ii,empty_indices)
                measurename = nii(jj).name;
                tractprofiles.(tractname).(measurename).profile = NaN(numnodes,1);
                tractprofiles.(tractname).(measurename).mean = NaN(1,1);
                tractprofiles.(tractname).(measurename).sd = NaN(1,1);
                tractprofiles.(tractname).x_coords = NaN(numnodes,1);
                tractprofiles.(tractname).y_coords = NaN(numnodes,1);
                tractprofiles.(tractname).z_coords = NaN(numnodes,1);
            else
                measurename = nii(jj).name;
                tractprofiles.(tractname).(measurename).profile = [];
                tractprofiles.(tractname).(measurename).mean = [];
                tractprofiles.(tractname).(measurename).sd = [];   
            end
        end
    end
end        

% Set up cell for csv
tract_profiles = cell(numnodes, length(nii));

for ifg = 1:length(fg_classified)
    fg_filename = strrep(strrep(fg_classified{ifg}.name,'.','_'), ' ', '');
    try
        if config.fiberbased == 0
            display 'volume based statistics'
            [SuperFiber, fgResampled] = dtiComputeSuperFiberRepresentation(fg_classified{ifg}, [], numnodes);
            for jj = 1:length(nii)
                if length(fg_classified{ifg}.fibers) < 6
                    display('too few streamlines. outputting profile of NaNs')
                    nii(jj).mean = NaN(numnodes,1);
                    nii(jj).std = NaN(numnodes,1);
                else
                    display(sprintf('computing %s',nii(jj).name));
                    [tract, ~, ~, ~, ~, ~, ~, ~, ~, ~, myValsFgSTD] = dtiComputeDiffusionPropertiesAlongFG_sd(fgResampled, nii(jj).data,[],[],numnodes,[],SuperFiber);
                    nii(jj).mean = tract;
                    nii(jj).std = myValsFgSTD;
                end
                clear tract myValsFgSTD;
            end
        else
            display 'fiber based statistics'
            [SuperFiber, fgResampled] = dtiComputeSuperFiberRepresentation(dtiXformFiberCoords(fg_classified{ifg}, inv(nii(jj).data.qto_xyz),'img'),[], numnodes); % convert fibergroup to the proper space
            for jj = 1:length(nii)
                if length(fg_classified{ifg}.fibers) < 6
                    display('too few streamlines. outputting profile of NaNs')
                    nii(jj).mean = NaN(numnodes,1);
                    nii(jj).std = NaN(numnodes,1);
                else
                    display(sprintf('computing %s',nii(jj).name));
                    tract = Compute_FA_AlongFG(fg, nii(jj).data, [], [], numnodes,SuperFiber);
                    nii(jj).mean = mean(tract,'omitnan');
                    nii(jj).std = std(tract,'omitnan');
                end
                clear tract myValsFgSTD;
            end
        end
        
        for jj = 1:length(nii)
            tract_profiles(:,jj,1) = num2cell(nii(jj).mean);
            tract_profiles(:,jj,2) = num2cell(nii(jj).std);
        end
        
        jj = 0;
        tps = {'mean','sd'};
        for jj = 1:length(nii)
            for tp = 1:length(tract_profiles(1,jj,:))
                T(:,(2*jj-2)+tp) = table(tract_profiles(:,jj,tp));
                T.Properties.VariableNames{(2*jj-2)+tp} = sprintf('%s_%s',char(nii(jj).name),tps{tp});
                T.Properties.VariableUnits{(2*jj-2)+tp} = sprintf('%s',nii(jj).units);
                %T.Properties.VariableUnits{jj+1} = sprintf('%s',nii(jj).units);
            end
        end
        
        % set information for superfiber coordinates for QA and informative
        % figures
        T.x_coords = SuperFiber.fibers{:}(1,:)';
        T.Properties.VariableUnits{length(T.Properties.VariableNames)} = 'mm';
        T.y_coords = SuperFiber.fibers{:}(2,:)';
        T.Properties.VariableUnits{length(T.Properties.VariableNames)} = 'mm';
        T.z_coords = SuperFiber.fibers{:}(3,:)';
        T.Properties.VariableUnits{length(T.Properties.VariableNames)} = 'mm';
        
        writetable(T, strcat('profiles/', fg_filename, '_profiles.csv'));
        
        for jj = 1:length(nii)
            % think of a better heuristic for this. right now, if a measure
            % ends with an 'e', it's just going to be skipped. reason this
            % works is because currently, the only measure I allow to be
            % run that ends in 'e' are the inverse measures (i.e.
            % fa_inverse).
            if ~strcmp(nii(jj).name(end),'e') 
                measurename = nii(jj).name;
                tractprofiles.(fg_filename).(measurename).profile = round(cell2mat(tract_profiles(:,jj,1)'),4);
                tractprofiles.(fg_filename).(measurename).mean = round(mean(tractprofiles.(fg_filename).(measurename).profile),4);
                tractprofiles.(fg_filename).(measurename).sd = round(std(tractprofiles.(fg_filename).(measurename).profile),4);
            end
        end
        
        % add coords to structure for product.json
        tractprofiles.(fg_filename).x_coords = SuperFiber.fibers{:}(1,:);
        tractprofiles.(fg_filename).y_coords = SuperFiber.fibers{:}(2,:);
        tractprofiles.(fg_filename).z_coords = SuperFiber.fibers{:}(3,:);

        if exist(config.ad)
            % AD
            analysisProfiles(nii(1).mean,fgResampled,nii(1).name,'Axial Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(1).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_ad.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('Axial Diffusivity');
            numfiles = numfiles + 1;
            % FA
            analysisProfiles(nii(3).mean,fgResampled,nii(3).name,'Fractional Anisotropy',[0.00, 1.00],[0 .25 .5 .75],numnodes,nii(2).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_fa.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('Fractional Anistropy');
            numfiles = numfiles + 1;
            % MD
            analysisProfiles(nii(5).mean,fgResampled,nii(5).name,'Mean Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(3).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_md.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('Mean Diffusivity');
            numfiles = numfiles + 1;
            % RD
            analysisProfiles(nii(7).mean,fgResampled,nii(7).name,'Radial Diffusivity',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(4).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_rd.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('Radial Diffusivity');
            numfiles = numfiles + 1;
        end
        
        if exist(config.ga)
            % ga
            analysisProfiles(nii(9).mean,fgResampled,nii(9).name,'Geodesic Anisotropy',[0.00, 1.00],[0 .25 0.5 0.75],numnodes,nii(5).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_ga.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('Geodesic Anisotropy');
            numfiles = numfiles + 1;
            % ak
            analysisProfiles(nii(11).mean,fgResampled,nii(11).name,'Axial Kurtosis',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(6).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_ak.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('Axial Kurtosis');
            numfiles = numfiles + 1;
            % mk
            analysisProfiles(nii(13).mean,fgResampled,nii(13).name,'Mean Kurtosis',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(7).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_mk.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('Mean Kurtosis');
            numfiles = numfiles + 1;
            % rk
            analysisProfiles(nii(15).mean,fgResampled,nii(15).name,'Radial Kurtosis',[0.00, 2.00],[0 .5 1 1.5],numnodes,nii(8).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_rk.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('Radial Kurtosis');
            numfiles = numfiles + 1;
        end

        if isfield(config,'ndi')
            % NDI
            analysisProfiles(nii(end_index-6+1).mean,fgResampled,nii(end_index-6+1).name,'NDI',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+1).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_ndi.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('NDI');
            numfiles = numfiles + 1;
            % ISOVF
            analysisProfiles(nii(end_index-6+3).mean,fgResampled,nii(end_index-6+3).name,'ISOVF',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+2).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_isovf.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('ISOVF');
            numfiles = numfiles + 1;
            % ODI
            analysisProfiles(nii(end_index-6+5).mean,fgResampled,nii(end_index-6+5).name,'ODI',[0 1.00],[0.25 .5 .75],numnodes,nii(end_index-6+3).units);
            json.images(numfiles).filename = strcat('images/',fg_filename,'_odi.png');
            json.images(numfiles).name = fg_filename;
            json.images(numfiles).desc = strcat('ODI');
            numfiles = numfiles + 1;
        end
        
    catch ME
        possible_error=1;
        failed_tracts = strcat(failed_tracts,fg_filename," ");
        save('profiles/error_messages.mat','ME');
    end
    
    if length(fg_classified{ifg}.fibers) < 6
        possible_error_lows=1;
        failed_tracts_lows = strcat(failed_tracts_lows,fg_filename," ");
        save('profiles/error_messages_lows.mat','failed_tracts_lows')
    end
    
    clf
    clear fgResampled SuperFiber;
end

fileID = fopen('numfiles.txt','w');
fprintf(fileID, '%d', numfiles-1); %matlab uses 1 based indexing
fclose(fileID);

message = struct;
if possible_error == 1
    message.type = 'error';
    message.msg = sprintf('ERROR: The following tracts have failed: %s',failed_tracts);
elseif possible_error_lows==1
    message.type = 'error';
    message.msg = sprintf('ERROR: The following tracts have too few streamlines: %s',failed_tracts_lows);
else
    message.type = 'success';
    message.msg = 'All tracts analysis profiles were created successfully';
end

product = struct;
product.brainlife = {message};
product.profiles = tractprofiles;

savejson('', product, 'product.json');
savejson('', json, fullfile('images.json'));

end