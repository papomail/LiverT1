%%CHECK DATA INTEGRITY & PARSE IT TO A PREPROCESS FOLDER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File Selection using Filter 
clear mFA_files filtered2
% mFA files
for patient=1:numel(pfolders)
    FAfiltered=FilterFiles(pfolders{patient},'FA');
    nfiltered=numel(FAfiltered);
    filtered=FAfiltered;
    if nfiltered>5
        count=1;
        for ii=1:nfiltered
            if ~isempty(strfind(FAfiltered(ii).name,'C.nii')) && isempty(strfind(FAfiltered(ii).name,char([filesep,'o'])))
                filtered2(count)=FAfiltered(ii);
                count=count+1;
            end
        end
        filtered=filtered2;
        nfiltered=numel(filtered);
    end
    
    
    if nfiltered<5
        clear Bfiles
        count=1;
        for ii=1:numel(FAfiltered);
            if ~isempty(strfind(FAfiltered(ii).name,'B.nii')) && isempty(strfind(FAfiltered(ii).name,char([filesep,'o'])))
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
                if ~exist(fullfile(outdir,[basename,'0001.nii.gz']),'file')
                    system(['cd ',outdir,'; fslsplit ',Bfiles(ii).name,' ',basename ]);
                end
                filtered(numel(filtered)+1)=subdir(fullfile(outdir,[basename,'0001.nii.gz']));
            end
        end
        
        
    end
    
    mFA_files{patient}=filtered;
    
    
    
   % check FA data is complete
   if numel(mFA_files{patient})==5
       disp(['mFA data of patient ',mat2str(patient),' seems to be complete.'])
   else
       error(['mFA data of patient ',mat2str(patient),' is incomplete... please check the dataset.'])
   end
   
   
end
 


%%  B1 files

clear B1_files
for patient=1:numel(pfolders)
    
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
    end
    
    B1_files{patient}=FilterFiles(pfolders{patient},'B1mapVol');
    disp(['B1map volume of patient ',mat2str(patient),' is complete.'])
end




%% Find FA (from protocol_name...)
clear FAarrays
for patient=1:numel(pfolders)
    for ii=1:5
        s=strfind(mFA_files{patient}(ii).name,'FA') ;
        ss=mFA_files{patient}(ii).name(s+2:s+3);
        sss=str2double(ss);
        if isnan(sss)
            ss=mFA_files{patient}(ii).name(s+2);
            sss=str2double(ss);
        end
        if sss==25
            sss=2.5;
        end
         FAarrays{patient}(ii)=sss;   
        
    end
end


%% Save file info as .mat (redundant)
save([fileparts(fileparts(pfolders{1})),filesep,'to_process_ALL.mat'],'FAarrays','B1_files','mFA_files');
for patient=1:numel(pfolders)
    FAarray=FAarrays{patient};
    B1_file=B1_files{patient};
    mFA_file=mFA_files{patient};
save([pfolders{patient},filesep,'to_process.mat'],'FAarray','B1_file','mFA_file');
end


%% Copy parsed files to a 'preprocessed' folder

for patient=1:numel(pfolders)
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
end