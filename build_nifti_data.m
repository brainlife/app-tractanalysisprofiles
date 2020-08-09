function nii = build_nifti_data(data,scale_value,value_units,inverse_units)
    nii = [];
    end_index = 0;
    for ii = 1:length(data)
        nii(end_index+1).name = char(extractBefore(data{ii}.name,strlength(data{ii}.name)-6));
        nii(end_index+1).data = niftiRead(fullfile(data{ii}.folder,data{ii}.name));
        nii(end_index+1).non_zero_index = find(nii(end_index+1).data.data(:,:,:) ~= 0);
        if strcmp(scale_value(ii),'true')
            if median(nii(end_index+1).data.data(nii(end_index+1).non_zero_index)) < 0.01
                nii(end_index+1).data.data = nii(end_index+1).data.data * 1000;
            end
        end
        nii(end_index+1).data_inv = 1./nii(end_index+1).data.data;
        nii(end_index+1).data_inv(~isfinite(nii(end_index+1).data_inv))=0;
        nii(end_index+1).units = value_units(ii)
        end_index = length(nii);
        nii(end_index+1).name = strcat(char(extractBefore(data{ii}.name,strlength(data{ii}.name)-6)),'_inverse');
        nii(end_index+1).data = nii(end_index).data;
        nii(end_index+1).data.data = nii(end_index).data_inv;
        nii(end_index+1).data_inv = nii(end_index).data_inv;
        nii(end_index+1).non_zero_index = nii(end_index).non_zero_index;
        nii(end_index+1).units = inverse_units(ii);
        end_index = length(nii);
    end
end