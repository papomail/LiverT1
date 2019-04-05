function out = ffiler(varargin)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% check if folders are already folderfiles



out.files=0;
out.skipmoco=0;
out.skippatient=0;
out.maskfile={};


if numel(varargin)<1
    message='Select the .nii files to do motion correction on, putting Reference image first and the Mask last. ';
    files=uipickfiles('FilterSpec',the_folder , 'REFilter','spirTSE|TR|0.nii\','Prompt',message,'Type', {'*.nii','NIFTI files'});
else
    options = struct(varargin{:});
    folder=options.folder;
    filter={options.filter};
    
    
    
    
    %% Check if NIFTI folder exists. If it doesnt convert dicom to nifti
    Folder_niftii = 'NIFTI';
    niftifolder=[folder,filesep,Folder_niftii];
    if ~exist(niftifolder,'dir')
        
        [~] = dcm2nii_recursive(folder,'',Folder_niftii); % convert the data to .nii
    end
    %%%%%%
    
    
    %% File selection
    for ii=1:numel(filter)
        myfilter=['*',filter{ii},'*'];
        files_st=subdir([folder,filesep,myfilter]);
        if ~isempty(files_st)
            files_temp{ii}={files_st.name};
        else
            files_temp{ii}={};
        end
    end
    
    %% Merge all files filtered
    files=[];
    for ii=1:numel(filter)
        files=cat(2,files,files_temp{ii});
    end
    
    %% Exclude some files  (T1MAP...)
    %%%%% dont choose T1MAPs
    string=strfind(files,'T1MAP');
    string=find(cellfun(@isempty,string));  %% choose files that are not T1MAPs
    files=files(string);
    
    % dont run in FLIRT is found
    string=strfind(files,'FLIRT');
    %string=find(~cellfun(@isempty,string));
    if any(~cellfun(@isempty,string))
        disp('FLIRT already done. Skipping MOCO...')
        string=find(cellfun(@isempty,string));  %% choose files that are not FLIRT
        out.skipmoco=1;
        %return
    end
    
    
    %Check if data complete
    if numel(files)~=5
        disp(['Data incomplete or inconsistent. Please review data in ',folder,'. Skipping Patient...']);
        out.skippatient=1;
        %return
    end
    
    %Check for liver mask
    
    if isfield(options,'maskfilter')
        maskfilter=options.maskfilter;
        %maskfile=dirrec(folder,[maskfilter,'*']);
        mymaskfilter=['*',maskfilter,'*'];
        maskfile_st=subdir([folder,filesep,mymaskfilter]);
        if ~isempty(maskfile_st)
            maskfile={maskfile_st.name};
        else
            maskfile={};
        end
        
        if isempty(maskfile)
            mymaskfilter=lower(mymaskfilter);
            maskfile_st=subdir([folder,filesep,mymaskfilter]);
            if ~isempty(maskfile_st)
                maskfile={maskfile_st.name};
            else
                maskfile={};
            end
            if isempty(maskfile)
                disp(['Can not find mask of liver boundary. Please review data in ',folder,'. Skipping MOCO...']);
                
            end
        end
        out.maskfile=maskfile;
    end
    % add mask flie to the list of files
    %files{numel(files)+1}=maskfile;
    
    
    
    out.files=files ;
    
    
end

