function out = parsedata(varargin)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if numel(varargin)<1
    message='Select the .nii files to do motion correction on, putting Reference image first and the Mask last. ';
    files=uipickfiles('FilterSpec',the_folder , 'REFilter','spirTSE|TR|0.nii\','Prompt',message,'Type', {'*.nii','NIFTI files'});
else
    options = struct(varargin{:});
    folder=options.folder;
    filter=options.filter;
    
    myfilter=[filter,'*'] ;
    files=dirrec([folder,filesep,'NIFTI'],myfilter);
    if isempty(files)
        myfilter=['*',myfilter];
        files=dirrec([folder,filesep,'NIFTI'],myfilter);
    end
    
    
    % dont choose T1MAPs
    string=strfind(files,'T1MAP');
    string=find(cellfun(@isempty,string));  %% choose files that are not T1MAPs
    files=files(string);
    
    % dont run in FLIRT is found
    string=strfind(files,'FLIRT');
    %string=find(~cellfun(@isempty,string));
    if any(~cellfun(@isempty,string))
        disp('FLIRT already done. Skipping MOCO...')
        out.skipmoco=1;
        return
    end
    
    
    %Check if data complete
    if numel(files)~=5
        disp(['Data incomplete or inconsistent. Please review data in ',folder,'. Skipping patient...']);
        out.skippatient=1;

        return
    end
    
    %Check for liver mask
    if isfield(options,'maskfilter')
        maskfilter=options.maskfilter;      
        maskfile=dirrec(folder,[maskfilter,'*']);
        if isempty(maskfile)
            maskfilter(1)=lower(maskfilter(1));
            maskfile=dirrec(folder,[maskfilter,'*']);
            if isempty(maskfile)
                disp(['Can not find mask of liver boundary. Please review data in ',folder,'. Skipping MOCO...']);
                        out.skipmoco=1;
                        return 

            end
        end
        out.maskfile=maskfile;
    end
    % add mask flie to the list of files
    %files{numel(files)+1}=maskfile;
end





end

