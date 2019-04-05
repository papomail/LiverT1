function VOL_MOCO_LIVER(varargin)

DO_SBS=0; %to perform slice by sclice moco (for IR data)

SKIP_4D_CHECK=1; % no need for this if preprocess folder contains the 3d-volumes (water).


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

use_schedule=0; %(By default do not use FLIRT schedule. It will be used if set in varargin )
if numel(varargin)<1
    message='Select the .nii files to do motion correction on, putting Reference image first and the Mask last. ';
    files=uipickfiles('FilterSpec',the_folder , 'REFilter','spirTSE|TR|0.nii\','Prompt',message,'Type', {'*.nii','NIFTI files'});
else
    options = struct(varargin{:});
    files={options.files};
    
    if isfield(options,'schedule')
        use_schedule=1;
        schedule=options.schedule;
    end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start MOCO
disp('VOLUME MOCO in progresss... please be patient');
mocotic=tic;



[mydir,~,~]=fileparts(files{1});
NumberLoaded=numel(files);

%split_ref=['cd ',dir,'; fslsplit ',refim,' aREF -z' ] ;

% Mask file
mask=files{NumberLoaded};

% image files
all_im=files(1:end-1);

%% STEP: choose one of the dimension (water) if mDixon data
clear waterfiles

if ~SKIP_4D_CHECK
    [~,dim4]=system(['fslval ',files{1},' dim4']);
    dim4=str2num(dim4);
    
    if dim4>1
        
        splittic=tic;
        for ii=1:numel(all_im)
            
            favol=files{ii};
            [~,b,ext]=fileparts(favol);
            if strcmp(ext,'.gz')
                waterfiles{ii}=['W_',b,'.gz'];
            else
                waterfiles{ii}=['W_',b,'.nii.gz'];
            end
            
            disp(['Spliting ', favol])
            cmd_slice_mov=['cd ',mydir,'; fslsplit ',favol, '; rm -f vol0000.nii.gz ; rm -f vol0002.nii.gz ; rm -f vol0003.nii.gz ; mv vol0001.nii.gz ',waterfiles{ii}] ;
            system(cmd_slice_mov);
            
        end
        splittime=num2str(round(toc(splittic)/60,1));
        disp(['It took ',splittime ,' minutes to split the volumes.',]);
    else
        waterfiles=  all_im;
    end
else
    waterfiles=  all_im;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Motion Correction (for mFA data)
[~,r,~]=fileparts(waterfiles{1});%% reference
refim1=waterfiles{1};
Ntomoco=numel(waterfiles)-1;
parfor ii=1:Ntomoco
    movim1=waterfiles{ii+1};
    [~,b,c]=fileparts(movim1);
        %outim1=[b,'_3DREG'];
        outim1=['3DREG_',b,c];
    moco3dtic=tic;
    disp(['Performing 3D moco of ', b, ' using ', r, ' as reference.'])
    if use_schedule
      cmd_3dmoco=['cd ',mydir,'; flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask,' -inweight ',mask,' -out ',outim1, ' -schedule ',schedule] ;
    else
      cmd_3dmoco=['cd ',mydir,'; flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask,' -inweight ',mask,' -out ',outim1] ;
    end
    %cmd_3dmoco=['flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask1,' -inweight ',mask1,' -out ',outim1,' -cost normmi -schedule /Users/patxi/Desktop/ytransonly.sch'] ;
    disp(cmd_3dmoco)
    system(cmd_3dmoco);
    moco3dtime=num2str(round(toc(moco3dtic)/60,1));
    disp([b ' done in ',moco3dtime ,' minutes. (',char(datetime),').']);
    %disp(datetime);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slice by Slice Motion Correction (for IR data)


if DO_SBS
%% STEP: Create the reference image.
switch refcase
    case 'fancy'
        
        dotempmean=['cd ',mydir,'; fsladd tempmean -m ',strjoin(all_im)] ; %create temporary mean-image
        system(dotempmean);
        
        weight=zeros(1:numel(all_im));
        for ii=1:numel(all_im)
            minusmean=['cd ',mydir,'; fslmaths tempmean -sub ',strjoin(all_im(ii)),' -abs minustempmean_',num2str(ii)] ; %substract mean-image to each image (and take the abs value)
            system(minusmean);
            
            stats=['cd ',mydir,'; fslstats -K ',mask,' minustempmean_',num2str(ii),' -m'];
            [~, meanval]=system(stats);
            weight(ii)=1/str2double(meanval);
            
            weighted=['cd ',mydir,'; fslmaths ',strjoin(all_im(ii)),' -mul ' ,num2str(weight(ii)),' weightedtemp',num2str(ii)] ; %substract mean-image to each image (and take the abs value)
            system(weighted);
        end
        
        [~, weightedim]=system(['ls ',mydir,filesep,'* |grep weightedtemp ']);
        goodmean=['cd ',mydir,'; fsladd myRef -m ',strjoin(strsplit(weightedim))] ; %create temporary mean-image
        system(goodmean);
        
        deletetemp=['cd ',mydir,'; find * | grep temp | xargs rm'] ; %% Be very carefull!!!
        system(deletetemp);
        
    case 'mean'
        domean=['cd ',mydir,'; fsladd myRef -m ',strjoin(all_im)] ; %create temporary mean-image
        system(domean);
        
    case 'max'
        dovol=['cd ',mydir,'; fslmerge -t allvolstemp ',strjoin(all_im)] ; %create volume
        system(dovol);
        domax=['cd ',mydir,'; fslmaths allvolstemp -Tmax myRef'] ; %max of the volume
        system(domax);
        deletetemp=['cd ',mydir,'; find * | grep temp | xargs rm'] ; %% Be very carefull!!!
        system(deletetemp);
end
%%


% Split into slices
if strcmp(refcase,'first')
    slice_ref=['cd ',mydir,'; fslslice ',files{1},' aREF' ] ;
else
    slice_ref=['cd ',mydir,'; fslslice myRef aREF' ] ;
end

system(slice_ref);
[~, refs]=system(['ls ',mydir,filesep,'* |grep aREF_slice |grep .gz']);
refs=strsplit(refs);
%refs=strcat([dir,filesep],refs);

slice_mask=['cd ',mydir,'; fslslice ',mask,' aMASK' ] ;
system(slice_mask);
[~, masks]=system(['ls ',mydir,filesep,'* |grep aMASK_slice |grep .gz']);
masks=strsplit(masks);
%masks=strcat([dir,filesep],masks);

%STEP: FLIRT of each slices
for ii=1:NumberLoaded-1
    movim_all=files{ii};
    
    cmd_slice_mov=['cd ',mydir,'; fslslice ',movim_all,' aMOV' ] ;
    system(cmd_slice_mov);
    
    [~, movs]=system(['ls ',mydir,filesep,'* |grep aMOV_slice |grep .gz |grep -v STEP']);
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
    
    [~, tomerge]=system(['ls ',mydir,filesep,'* |grep aMOV_slice |grep .gz |grep STEP']);
    tomerge=strsplit(tomerge);
    %tomerge=string(tomerge);
    tomerge=strjoin(tomerge);
    
    system(['fslmerge -a ',movim_all,'_FLIRT ',char(tomerge)]);
end

end % endof SBS moco

moco_time=num2str(round(toc(mocotic)/60,1));
disp(['Time for MOCO was ',moco_time,' minutes']);
end

