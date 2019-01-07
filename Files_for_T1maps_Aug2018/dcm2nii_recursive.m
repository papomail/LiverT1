function [path_nii] = dcm2nii_recursive(dirnames,source,destination)


    for ik=1:numel(dirnames)
        pathname=dirnames{ik};

             pathname_slash = [pathname, slsh];
             path_nii{ik} = [pathname_slash,''];

                mkdir([path_nii{ik},destination])
    %        convert from DCM to NII (if dcm are multi-frame, it will give single-frame .nii) 
             dicm2nii([path_nii{ik},source], [path_nii{ik},destination], 'nii'); % 
    end
end

