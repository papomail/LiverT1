function [ files] = MOCO_LIVER_Jan2019(varargin)

%refcase='mean'
refcase='first';
%refcase='fancy'
%refcase='max'


%set FSL invironment
setenv('FSLDIR','/usr/local/fsl/')
setenv('FSLOUTPUTTYPE','NIFTI_GZ')


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[the_folder] = pathsetting('Desktop');

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
        return
    end
    
    
    %Check if data complete
    if numel(files)~=5
        disp(['Data incomplete or inconsistent. Please review data in ',folder,'. Skipping MOCO...']);
        return
    end
    
    %Check for liver mask
    maskfile=dirrec(folder,'Seg*');
    if isempty(maskfile)
        maskfile=dirrec(folder,'seg*');
        if isempty(maskfile)
        disp(['Can not find mask of liver boundary. Please review data in ',folder,'. Skipping MOCO...']);
        return
        end
    end
    
    
    % add mask flie to the list of files
    files{numel(files)+1}=maskfile;
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start MOCO
disp('MOCO in progresss... please be patient');
mocotic=tic;



[dir,~,~]=fileparts(files{1});
NumberLoaded=numel(files);

%split_ref=['cd ',dir,'; fslsplit ',refim,' aREF -z' ] ;

% Mask file
mask=files{NumberLoaded};

% image files
all_im=files(1:end-1);

%% STEP: choose one of the dimension (water) if mDixon data


%% STEP: Create the reference image.
switch refcase
    case 'fancy'
        
        dotempmean=['cd ',dir,'; fsladd tempmean -m ',strjoin(all_im)] ; %create temporary mean-image
        system(dotempmean);
        
        weight=zeros(1:numel(all_im));
        for ii=1:numel(all_im)
            minusmean=['cd ',dir,'; fslmaths tempmean -sub ',strjoin(all_im(ii)),' -abs minustempmean_',num2str(ii)] ; %substract mean-image to each image (and take the abs value)
            system(minusmean);
            
            stats=['cd ',dir,'; fslstats -K ',mask,' minustempmean_',num2str(ii),' -m'];
            [~, meanval]=system(stats);
            weight(ii)=1/str2double(meanval);
            
            weighted=['cd ',dir,'; fslmaths ',strjoin(all_im(ii)),' -mul ' ,num2str(weight(ii)),' weightedtemp',num2str(ii)] ; %substract mean-image to each image (and take the abs value)
            system(weighted);
        end
        
        [~, weightedim]=system(['ls ',dir,filesep,'* |grep weightedtemp ']);
        goodmean=['cd ',dir,'; fsladd myRef -m ',strjoin(strsplit(weightedim))] ; %create temporary mean-image
        system(goodmean);
        
        deletetemp=['cd ',dir,'; find * | grep temp | xargs rm'] ; %% Be very carefull!!!
        system(deletetemp);
        
    case 'mean'
        domean=['cd ',dir,'; fsladd myRef -m ',strjoin(all_im)] ; %create temporary mean-image
        system(domean);
        
    case 'max'
        dovol=['cd ',dir,'; fslmerge -t allvolstemp ',strjoin(all_im)] ; %create volume
        system(dovol);
        domax=['cd ',dir,'; fslmaths allvolstemp -Tmax myRef'] ; %max of the volume
        system(domax);
        deletetemp=['cd ',dir,'; find * | grep temp | xargs rm'] ; %% Be very carefull!!!
        system(deletetemp);
end
%%


% Split into slices
if strcmp(refcase,'first')
    slice_ref=['cd ',dir,'; fslslice ',files{1},' aREF' ] ;
else
    slice_ref=['cd ',dir,'; fslslice myRef aREF' ] ;
end

system(slice_ref);
[~, refs]=system(['ls ',dir,filesep,'* |grep aREF_slice |grep .gz']);
refs=strsplit(refs);
%refs=strcat([dir,filesep],refs);

slice_mask=['cd ',dir,'; fslslice ',mask,' aMASK' ] ;
system(slice_mask);
[~, masks]=system(['ls ',dir,filesep,'* |grep aMASK_slice |grep .gz']);
masks=strsplit(masks);
%masks=strcat([dir,filesep],masks);

%STEP: FLIRT of each slices
for ii=1:NumberLoaded-1
    movim_all=files{ii};
    
    cmd_slice_mov=['cd ',dir,'; fslslice ',movim_all,' aMOV' ] ;
    system(cmd_slice_mov);
    
    [~, movs]=system(['ls ',dir,filesep,'* |grep aMOV_slice |grep .gz |grep -v STEP']);
    movs=strsplit(movs);
    %movs=strcat([dir,filesep],movs);
    
    %remove_temp_movs='rm
    
    for jj=1:numel(movs)-1
        mask1=masks{jj};
        
        outim1=[movs{jj},'_aSTEP1'];
        refim1=refs{jj};
        movim1=movs{jj};
        
        disp(['Registering volume',num2str(ii),' slice',num2str(jj)])
        if strcmp(refcase,'first') && ii==1
            system(['fslmaths ',refim1,' -mul 1 ',outim1]); %so that it does not registre to itself
        else
            %command=['flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask1,' -out ',outim1,' -2D -cost normmi '] ;
            %command=['flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask1,' -out ',outim1,' -2D -cost normmi  -searchrx -5 5 -searchry -5 5 -searchrz -5 5'] ;
            command=['flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask1,' -inweight ',mask1,' -out ',outim1,' -cost normmi -schedule /Users/patxi/Desktop/ytransonly.sch'] ;
            
            system(command);
        end
        
    end
    
    [~, tomerge]=system(['ls ',dir,filesep,'* |grep aMOV_slice |grep .gz |grep STEP']);
    tomerge=strsplit(tomerge);
    %tomerge=string(tomerge);
    tomerge=strjoin(tomerge);
    
    system(['fslmerge -a ',movim_all,'_FLIRT ',char(tomerge)]);
end



moco_time=num2str(round(toc(mocotic)/60,1));
disp(['Time for MOCO was ',moco_time,' minutes']);
end

