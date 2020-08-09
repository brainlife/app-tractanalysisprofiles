function nii = build_nifti_data(data,scale_value,value_units,inverse_units)
    for ii = 1:length(data)
        nii(ii).name = char(extractBefore(data(ii).name,strlength(data(ii).name)-6));
        nii(ii).data = niftiRead(fullfile(data(ii).folder,data(ii).name));
        nii(ii).non_zero_index = find(nii(ii).data.data(:,:,:) ~= 0);
        if strcmp(scale_value(ii),'true')
            if median(nii(ii).data.data(nii(ii).non_zero_index)) < 0.01
                nii(ii).data.data = nii(ii).data.data * 1000;
            end
        end
        nii(ii).data_inv = 1./nii(ii).data.data;
        nii(ii).data_inv(~isfinite(nii(ii).data_inv))=0;
        nii(ii).units = value_units(ii)
        end_index = length(nii);
        nii(end_index+ii).name = strcat(char(extractBefore(data(ii).name,strlength(data(ii).name)-6)),'_inverse');
        nii(end_index+ii).data = nii(ii).data;
        nii(end_index+ii).data.data = nii(ii).data_inv;
        nii(end_index+ii).data_inv = nii(ii).data_inv;
        nii(end_index+ii).units = inverse_units(ii);
    end
end