%%CHECK DATA INTEGRITY & PARSE IT TO A PREPROCESS FOLDER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File Selection using Filter 
clear mFA_files filtered2
% mFA files
NIFTIfolder=[pfolders{patient},filesep,'NIFTI'];
    FAfiltered=FilterFiles(NIFTIfolder,'FA');
    nfiltered=numel(FAfiltered);
    filtered=FAfiltered;
    
    
    if nfiltered>5
        count=1;
        for ii=1:nfiltered
            if ~isempty(strfind(FAfiltered(ii).name,'0003.nii.gz'))
                filtered2(count)=FAfiltered(ii);
                count=count+1;
            end
        end
        if exist('filtered2','var')
        filtered=filtered2;
        nfiltered=numel(filtered);
        end
    end
    
    
    
    
    
    if nfiltered>5
        count=1;
        for ii=1:nfiltered
            if ~isempty(strfind(FAfiltered(ii).name,'C.nii')) && isempty(strfind(FAfiltered(ii).name,char([filesep,'o']))) && isempty(strfind(FAfiltered(ii).name,'Split')) && isempty(strfind(FAfiltered(ii).name,'000'))
                filtered2(count)=FAfiltered(ii);
                count=count+1;
            end
        end
         if exist('filtered2','var')
        filtered=filtered2;
        nfiltered=numel(filtered);
         end
    end
    
    
    if nfiltered<5
        clear Bfiles
        count=1;
        for ii=1:numel(FAfiltered);
            if ~isempty(strfind(FAfiltered(ii).name,'B.nii')) && isempty(strfind(FAfiltered(ii).name,char([filesep,'o'])))  && isempty(strfind(FAfiltered(ii).name,'Split')) && isempty(strfind(FAfiltered(ii).name,'000'))
                Bfiles(count)=FAfiltered(ii);
                count=count+1;
            end
        end
        for ii=1:numel(Bfiles)
            [~, bb]=system(['fslinfo ',Bfiles(ii).name]);
            dim4=str2double(bb(85:95));
            if dim4>1
                [outdir,basename,~]=fileparts(Bfiles(ii).name);
                if strfind(basename,'.nii')
                    [~,basename,~]=fileparts(basename);
                end
                basename=char([basename,'Split']);
                if ~exist(fullfile(outdir,[basename,'0001.nii.gz']),'file')
                    system(['cd ',outdir,'; fslsplit ',Bfiles(ii).name,' ',basename ]);
                end
                filtered(numel(filtered)+1)=subdir(fullfile(outdir,[basename,'0001.nii.gz']));
            end
        end 
        
    end
    
    
    
    mFA_files{patient}=filtered;
    
    
    
    %% check FA data is complete and split volumes if necessary
    if numel(mFA_files{patient})==5
        disp(['mFA data of patient ',mat2str(patient),' seems to be complete.'])
        clear singleVolFa
        parfor ii=1:5
            
            temp_fa=mFA_files{patient}(ii).name;
            [~, bb]=system(['fslinfo ',temp_fa]);
            dim4=str2double(bb(85:95));
            if dim4>1
                [outdir,basename,~]=fileparts(temp_fa);
                if strfind(basename,'.nii')
                    [~,basename,~]=fileparts(basename);
                end
                basename=char([basename,'Split']);
                
                if ~exist(fullfile(outdir,[basename,'0001.nii.gz']),'file')
                    system(['cd ',outdir,'; fslsplit ',temp_fa,' ',basename ]);
                end
                disp(['Splitting 4D volume ',num2str(ii),' of 5']);
                singleVolFa(ii)=subdir(fullfile(outdir,[basename,'0003.nii.gz']));
            else
                singleVolFa(ii)=mFA_files{patient}(ii);
            end
            
        end
        mFA_files{patient}=singleVolFa;
        
    else
        error(['mFA data of patient ',mat2str(patient),' is incomplete... please check the dataset.'])
    end
    
    
    
    


%%  B1 files

clear B1_files

    myb1file=FilterFiles(NIFTIfolder,'B1mapVol');
    if isempty(FilterFiles(pfolders{patient},'B1mapVol'))
    b1slices=FilterFiles(pfolders{patient},'B1');
        if numel(b1slices)==20
            outdir=fileparts(b1slices(1).name);
            b1list=[];
            for ii=1:15,
                b1list=cat(2,b1list,[b1slices(ii).name,' ']);
            end
            system(['cd ',outdir,'; fslmerge -z B1mapVol ',b1list ]);
        else
            error(['B1 map of patient ',mat2str(patient),' seems to be incomplete... please check the dataset.'])
        end
        myb1file=FilterFiles(NIFTIfolder,'B1mapVol');
    end
    
    B1_files{patient}=myb1file;
    disp(['B1map volume of patient ',mat2str(patient),' is complete.'])








%% Save file info as .mat (redundant)
save([fileparts(fileparts(pfolders{1})),filesep,'to_process_ALL.mat'],'B1_files','mFA_files');

    %FAarray=FAarrays{patient};
    B1_file=B1_files{patient};
    mFA_file=mFA_files{patient};
save([pfolders{patient},filesep,'to_process.mat'],'B1_file','mFA_file');



%% Copy parsed files to a 'preprocessed' folder


    [a,b,~]=fileparts(pfolders{patient});
    [~,bb,~]=fileparts(a);
    prepfolder=fullfile(pfolders{patient},'preprocessed');
    mkdir(prepfolder);
    B1_file=B1_files{patient};
    disp(['Dataset ',char([bb,filesep,b]),' (',num2str(patient),' of ',num2str(numel(pfolders)),'):']);
    disp('Copying B1map to preprocessed folder')
    copyfile(B1_file.name, prepfolder);
    for ii=1:5
        mFA_file=mFA_files{patient}(ii);
        disp(['Copying ',num2str(ii),' of 5 mFA files to preprocessed folder'])
        copyfile(mFA_file.name, prepfolder);
    end
