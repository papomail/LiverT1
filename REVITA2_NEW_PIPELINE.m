%% T1 maps for REVITA2 study
%set FSL invironment
setenv('FSLDIR','/usr/local/fsl/')
setenv('FSLOUTPUTTYPE','NIFTI_GZ')

fclose all
format long
clear all

%pathtoadd=genpath(pwd);
%addpath(pathtoadd)
%% Options
BATCH_PROCCESS=1

SAVE_NII=1

SBS_MOCO=0 %% To apply slice by slice motion correction (useful for IR data)


DO_MATLAB_3D_MOCO=0 % To apply extra motion correction as a rigid volume with matlab (useful...?)

                    DO_SMOOTH_B1=1
                    EXTREME_SMOOTH_B1=1
                    DO_NOT_USE_B1_DATA=0  %%!!!

MERGEZSLICES=1 %% merges zslices in pairs to obtain similar slice thickness as in the IR_TSE sequence


fit_type='multi_fa'
use_b1_in_fit=1


%% Set original path
[the_folder] = pathsetting('Desktop/REVITA22');

%%  Patient Selection
message = 'Select patient folders. (One folder per visit)';
pfolders=uipickfiles('FilterSpec',the_folder,'Prompt',message);
%load('/Users/patxi/Desktop/REVITA22/pfolders.mat')


clear patientfolder patientvisit
[a,~,~]=fileparts(pfolders{1});
studyfolder=fileparts(a);


batchtic=tic;


processed=0; 

 %% log
    logfile=[studyfolder,filesep,'log_',mydate,'.txt'];
    diary(logfile)
    
%% Parse data
%parseout=parsedata('folder',pfolders{patient},'filter',{'W*FA10.nii.gz','3DREG'},'maskfilter','Seg');  %%  REVITA2 3DREGISTRED filenames


%% NIFTI conversion
parfor patient=1:numel(pfolders)
    
    if ~exist(fullfile(pfolders{patient},'NIFTI'),'dir')
    
    patienttic=tic;
    disp(' ');
    disp(' ');
    disp(' ');
    disp('**********************************************************************');
    disp(['Starting NIFTI conversion of ',pfolders{patient},':']);
    disp(' ');
    subjectFolder=pfolders{patient};
    NIFTIFolder = fullfile(pfolders{patient},'NIFTI');

    mkdir(NIFTIFolder);
    %unix([ '/Applications/osx/dcm2nii -b dcm2nii_rescale.ini -a N -g Y -x N -o ' NIFTIFolder , ' ' subjectFolder]);
     %system([ '/Applications/osx/dcm2nii -b ',studyfolder,filesep,'dcm2nii_rescale.ini -a Y -g Y -x N -o ' NIFTIFolder , ' ' subjectFolder]);
     %system([ '/Applications/dcm2niix_25-Nov-2018_mac/dcm2niix -p Y -o ' NIFTIFolder , ' ' subjectFolder]);
         system([ '/Applications/osx/dcm2nii -b /Users/patxi/Desktop/REVITA22/dcm2nii_rescale.ini -a Y -g Y -x N -o ' NIFTIFolder , ' ' subjectFolder]);

    end

end



 %% Load DICOMs (without NIFTI conversion)
% warning off
% dinfo = dparse_folders(pfolders{patient}, 'FA' ) ;
% warning on
% [vol, matp] = d2mat(dinfo,{'slice','fa','wfio'},'wfio',[1],'op','fp') ;
% svol=size(vol);
% 
% warning off
% dinfoB1 = dparse_folders('/Users/patxi/Desktop/REVITA22/05-006/Visit2','B1') ;
% warning on
% [volB1, matpb1] = d2mat(dinfoB1,{'slice','itype'},'itype',[10],'op','dv') ;
% svolB1=size(volB1);
% 
% 
% 
% 
% figure;
% sl=round(svol(3)/2);
% imagesc(reshape(vol(:,:,sl,:),[svol(1),svol(2)*svol(4)]))
% title('mFA water data')
% 
% figure;
% slb1=round(svolB1(3)/2);
% imagesc(volB1(:,:,slb1,:))
% title('B1 data')
% 

MYARRAY=1:numel(pfolders);
%MYARRAY=[1:22,24:48]
%% Data selection
clear B1parsed B1maps Segmentation
for patient=MYARRAY%1:numel(pfolders)
    prepfolder=fullfile(pfolders{patient},'preprocessed');
    if ~exist(prepfolder,'dir')
    PARSE_TO_PREPROCESS_SINGLE
    end
    FAparsed{patient}=FilterFiles(prepfolder,'FA');
    B1parsed{patient}=FilterFiles(prepfolder,'B1');
    Segmentation{patient}=FilterFiles(fileparts(prepfolder),'Seg');
    if isempty(Segmentation{patient})
        Segmentation{patient}=FilterFiles(fileparts(prepfolder),'seg');
    end
end

%% Segmentation to match mFA-data orientation
clear Segmentation_OK
for patient=MYARRAY%1:numel(pfolders)
    Segmentation_OK{patient}=FilterFiles(pfolders{patient},'OK_');
    if isempty( Segmentation_OK{patient})
        template=FAparsed{patient}(1).name;
        source=Segmentation{patient}.name;
        [a,b,c]=fileparts(source);
        result=fullfile(a,['OK_',b,c]);
        nii_xform(source, template, result)
        disp(['Matching mask orientation to mFA-data. [',num2str(patient),' of ',num2str(numel(pfolders)),']'])
         Segmentation_OK{patient}=FilterFiles(pfolders{patient},'OK_');
    end
    disp(['Mask orientation OK. [',num2str(patient),' of ',num2str(numel(pfolders)),']'])
end


%% FA data sorting
%% Find FA (from protocol_name...)
clear FAarrays FAsorted
for patient=MYARRAY%1:numel(pfolders)
    for ii=1:5
        s=strfind(FAparsed{patient}(ii).name,'FA') ;
        ss=FAparsed{patient}(ii).name(s+2:s+3);
        sss=str2double(ss);
        if isnan(sss)
            ss=FAparsed{patient}(ii).name(s+2);
            sss=str2double(ss);
        end
        if sss==25
            sss=2.5;
        end
         FAarrays{patient}(ii)=sss; 
    end
    [~, per]=sort(FAarrays{patient});
    FAsorted{patient}=FAparsed{patient}(per);
end


%% MoCo   (do not use 'parfor' here, VOL_MOCO_LIVER will call 'parfor')
clear FAmocoed_temp FAmocoed
for patient=MYARRAY%1:numel(pfolders)
    prepfolder=fullfile(pfolders{patient},'preprocessed');

    FAmocoed_temp=FilterFiles(prepfolder,'3DREG');
     if isempty(FAmocoed_temp)
        filesnmask={FAsorted{patient}.name};
        mask={Segmentation_OK{patient}.name};
        filesnmask(numel(filesnmask)+1)=mask;
        
        %myschedule='/Users/patxi/Desktop/REVITA22/flirtsch/xyztrans.sch';
        myschedule='/Users/patxi/Desktop/REVITA22/flirtsch/yztrans.sch';
        %myschedule='/Users/patxi/Desktop/REVITA22/flirtsch/ytransonly.sch';
        disp(' ')
        disp(' ')
        disp(['Starting MoCo for ',pfolders{patient},' [',num2str(patient),' of ',num2str(numel(pfolders)),']'])

        VOL_MOCO_LIVER('files',filesnmask,'schedule',myschedule);
        FAmocoed_temp=FilterFiles(prepfolder,'3DREG');
    end
    
     disp(['MoCoed data identified for ',pfolders{patient},' [',num2str(patient),' of ',num2str(numel(pfolders)),']'])
     
         FAmocoed{patient}(1)=FilterFiles(prepfolder,'FA25');
            FAmocoed{patient}(2:5)=FilterFiles(prepfolder,'3DREG');
     
     
end



%% FA sorting of the MOCOed data

clear SortedFAarrays FA_OK
%parfor patient=1:numel(pfolders)
for patient=MYARRAY%1:numel(pfolders)

    for ii=1:5
        s=strfind(FAmocoed{patient}(ii).name,'FA') ;
        ss=FAmocoed{patient}(ii).name(s+2:s+3);
        sss=str2double(ss);
        if isnan(sss)
            ss=FAmocoed{patient}(ii).name(s+2);
            sss=str2double(ss);
        end
        if sss==25
            sss=2.5;
        end
         FAarrays{patient}(ii)=sss; 
    end
    [SortedFAarrays{patient}, per]=sort(FAarrays{patient});
    FA_OK{patient}=FAmocoed{patient}(per);
    
    disp(['Data sorted and ready to process in ',pfolders{patient},' [',num2str(patient),' of ',num2str(numel(pfolders)),']'])

end




%% Load data in MATLAB for T1 processing
% For each patient use data in:
% FA_OK{patient} for the mFA data
% B1parsed{patient} for the B1 maps
% Segmentation_OK{patient} for the mask (used for MoCo)
batchtic=tic;

%for patient=1:numel(pfolders)   %% (do not use 'parfor' here, multiFAfit will call 'parfor')
for patient=MYARRAY%1:numel(pfolders)
  
    disp(' ');
    disp(' ');
    disp(' ');
    disp('**********************************************************************');
    
    
    %% check if data already proccessed
    if  DO_SMOOTH_B1==0
        checkfor='_mFAfit_T1MAPS.nii.gz';
    end
  
    if  DO_SMOOTH_B1 && EXTREME_SMOOTH_B1
        checkfor='_mFAfit_T1MAPS_EXTREMEsmoothB1.nii.gz';
    end
    
    if  DO_SMOOTH_B1 && ~EXTREME_SMOOTH_B1
        checkfor='_mFAfit_T1MAPS_smoothB1.nii.gz';
    end
    
    mycheck=FilterFiles(pfolders{patient},checkfor);
    
    if ~isempty(mycheck)
        disp([ 'Data already processed with current paramenters. Skipping ',pfolders{patient},' [',num2str(patient),' of ',num2str(numel(pfolders)),']'])
        continue
    end
    
   %% 
    patienttic=tic;
    
    
    
    
    %% log
    logfile=[studyfolder,filesep,'log_',mydate,'.txt'];
    diary(logfile)
    
    % Load mFA data
    prepfolder=fullfile(pfolders{patient},'preprocessed');
    num_dataSets=numel(FA_OK{patient});
    clear dataT1_init niiSS niiRI
    disp([ 'Loading data of ',pfolders{patient},' [',num2str(patient),' of ',num2str(numel(pfolders)),']'])
    for ii = 1:num_dataSets  %% ii loops over the 5 mFA datasets per patient
        data_cell=nii_tool('load',FA_OK{patient}(ii).name);
        datafromnifti=double(data_cell.img);
        %%%% Scaling .nii for matlab to display float point
        niiSS(ii) = double(data_cell.hdr.scl_slope);
        niiRI(ii)= double(data_cell.hdr.scl_inter);
        dataT1_init(:,:,:,ii)=  datafromnifti*niiSS(ii)+niiRI(ii);
    end
    
    % Load B1 data
    B1map_data=nii_tool('load',B1parsed{patient}.name);
    if size(B1map_data.img,4)>1
        B1map_data.img(:,:,:,2:6)=[];
        B1map_data=nii_tool('update',B1map_data);
    end
    disp( 'Interpolating B1maps to match mFA data...')
    B1map_data_int=nii_xform(B1map_data,data_cell);
    B1map_data_int=B1map_data_int.img;
    %%%% Scaling .nii for matlab to display float point
    niiSS_b0 = double(B1map_data.hdr.scl_slope);
    niiRI_b0= double(B1map_data.hdr.scl_inter);
    B1map_data_init=double(B1map_data_int*niiSS_b0+niiRI_b0);
    
    
    % Load Segmentation file (for broad ROI selection)
    mask_file=nii_tool('load',Segmentation_OK{patient}.name);
    
    %% Rotate images (matlab opens them 90Clockwise)
    reorientdataT1=rot90(dataT1_init); %% rotate image 90Anti-clockwise
    B1map_data_init=rot90(B1map_data_init); %% rotate image 90A
    mask=rot90(mask_file.img);
    
    fs1=size(reorientdataT1,1);
    fs2=size(reorientdataT1,2);
    fs3=size(reorientdataT1,3);
    
    
     %% In-plane ROI selection for fitting
        
        if BATCH_PROCCESS
            if ~isempty(mask)
                %mym=nii_tool('load',char(mask));
                %mm=rot90(permute(mym.img,[1 3 2]));
                smm=nanmean(mask,3);
                mx=find(nanmean(smm,1));
                mx1=mx(1);
                my=find(nanmean(smm,2));
                my1=my(1);
                position=[round(mx1-0.1*fs1)   round(my1-0.1*fs1)   round(0.75*fs1)   round(0.75*fs2)];
                                
               % position=[50   50   200   200]; %for testing

                if position(1)<1
                    position(1)=1;
                end
                if position(2)<1
                    position(2)=1;
                end
                if position(3)>fs1-position(1)
                    position(3)=fs1-position(1);
                end
                if position(4)>fs2-position(2)
                    position(4)=fs2-position(2);
                end
                
            else
                position=[round(0.15*fs1)   round(0.23*fs2)   round(0.4*fs1)   round(0.5*fs2)]
            end
        else
            figure, imagesc(reorientdataT1(:,:,round(size(reorientdataT1,3)/2),1,1))
            title('Please draw a rectangle in the liver where fitting is required. [double clic to finish]')
            h = imrect;
            position = round(wait(h));
            close
        end
        
        xrange=position(1):position(1)+position(3);
        yrange=position(2):position(2)+position(4);
    
    %% Slice selection for fitting
    % [1/8 to 7/8] of the FOV in the Z direction
    slrange=[round(fs3/8):round(7*fs3/8)];
    if mod(numel(slrange),2)
        slrange=slrange(1:end-1);
    end
    
    %slrange=[60:61] %for testing
    
    dataT1=reorientdataT1(yrange,xrange,slrange,:);
    dataB1=B1map_data_init(yrange,xrange,slrange);
    
    
    if MERGEZSLICES
        disp('Merging zslices in pairs to obtain similar slice thickness as in the IR_TSE sequence')
        dataT1odd=dataT1(:,:,1:2:end,:);%% merge zslices in pairs
        dataT1even=dataT1(:,:,2:2:end,:);
        dataT1merg=(dataT1odd+dataT1even)/2;
        dataT1=dataT1merg;
        
        dataB1odd=dataB1(:,:,1:2:end);%% merge zslices in pairs
        dataB1even=dataB1(:,:,2:2:end);
        dataB1merg=(dataB1odd+dataB1even)/2;
        dataB1=dataB1merg;
        clear  dataT1_mer dataB1merg dataB1even dataB1odd dataT1even dataT1odd
    end
    
    
 
    
    %figure,imagesc(dataT1(:,:,22,1))
      %  figure,imagesc(dataB1(:,:,22))
       % figure,imagesc(reorientdataT1(:,:,60,1))
    %% B1 Smoothing
    if DO_SMOOTH_B1
        if EXTREME_SMOOTH_B1
            disp('Applying EXTREME smoothing to B1 field map.')
            dataB1=smooth3(dataB1,'box',[51,51,21]);
        else
            disp('Applying smoothing to B1 field map.')
            dataB1=smooth3(dataB1,'box',[21,21,11]);
        end
    end
    
    %% T1 Fitting 
    disp(' ')
    disp('Starting Multi FA fit')
    tic
    TR=4;
    [fitparam] = multiFAfit(dataT1,SortedFAarrays{patient},TR,dataB1); %fit with B1
    toc
    xvals=SortedFAarrays{patient};
    
    
     %% R^2 calculation
        myfval=fitparam.fval;
        [ R2adjusted ] = R2calc( dataT1,xvals,myfval );
        % mask out data that fit not a good fit.
        do_mask_based_on_R2=0
        if do_mask_based_on_R2
            myr2lim=0.80;
            rsquarefilt=R2adjusted>myr2lim;
            rsquarefilt=double(rsquarefilt);
            rsquarefilt(rsquarefilt==0)=nan;
        else
            rsquarefilt=1;
            myr2lim=0;
        end
        
        T1s=fitparam.RelaxTime.*rsquarefilt;
        M0s=fitparam.M0.*rsquarefilt;
        %myfval=myfval.*rsquarefilt;
        

       
        % reshape for plots
        T1maps=reshape(T1s,size(T1s,1),size(T1s,2)*size(T1s,3));
        M0maps=reshape(M0s,size(T1s,1),size(T1s,2)*size(T1s,3));

        if ~BATCH_PROCCESS
            % plot T1maps
            caxmax=2;
            caxmin=0;
            
            figure,
            imagesc((T1maps/1000)),title({'T1 map [s]',['(Masked with Adjusted R^2 > ',num2str(myr2lim),')']}), colorbar
            colormap(jet)
            caxis([caxmin caxmax])
        end
     
        
        %% SAVE .nii MAPS
        if SAVE_NII==1
            
           
            T1MAPS_folder=fullfile(pfolders{patient},'T1_MAPS');
            mkdir(T1MAPS_folder);
            
            my_filename=FA_OK{patient}(1).name
            [~, the_filename,ext] = fileparts(my_filename);
            if strcmp(ext,'.gz')
                [~, the_filename,ext2] = fileparts(the_filename);
            end
            
          
            
            
            tosave_T1maps=T1s;
            tosave_T1maps(tosave_T1maps>5000)=0;
            tosave_T1maps(tosave_T1maps<1)=0;
            tosave_T1maps=int16(tosave_T1maps);
            
            
            s1 = size(tosave_T1maps,1);
            s2 = size(tosave_T1maps,2);
            s3 = size(tosave_T1maps,3);
            s4 = size(tosave_T1maps,4);
            
            %load header from the original T1data & create the T1maps with the same header
            %niiheader = load_untouch_nii(filepath_T1data); %load a cest scan
            mynii=data_cell;
            
            %get the rid of scaling slope (if any)
            mynii.hdr.scl_slope=1;
            mynii.hdr.scl_inter=0;
            
            switch fit_type
            
                case 'multi_fa'
                    
                    fullsizeim=double(rot90((data_cell.img))); %open original image to substitute with t1maps values (xray effect)
                    if size(fullsizeim,4)>1
                        fullsizeim=fullsizeim(:,:,:,2);%take the water only image
                    end
                    myint=fullsizeim(:,:,round(size(fullsizeim,3)/2),1,1);
                    myint(myint==0)=nan;
                    fullsizeim=fullsizeim*300/nanmedian(myint(:)); %fill the original image with t1maps values (xray effect)
                    
                    
                    if MERGEZSLICES
                        tosave_T1maps2=zeros(size(tosave_T1maps,1),size(tosave_T1maps,2),size(tosave_T1maps,3)*2); %put every fittes slice into 2 original-thinkness slices
                        tosave_T1maps2(:,:,1:2:end)=tosave_T1maps;
                        tosave_T1maps2(:,:,2:2:end)=tosave_T1maps;
                        
                        fullsizeim(yrange,xrange,slrange)=tosave_T1maps2;
                    else
                        fullsizeim(yrange,xrange,slrange)=tosave_T1maps;
                    end
                    
                    fullsizeim=rot90(fullsizeim,-1); %Re-orient to initial orientation
                    mynii.img=fullsizeim;
                    
                    
                    
                    %make sure to go one dimension less as we now have maps
                    %mynii.hdr.dime.dim(2:5) = [size(fullsizeim),1];
                    %mynii.hdr.dime.pixdim(5) = 0;
                    mynii=nii_tool('update',mynii);
                    %make sure output is float
                    %niiheader.hdr.dime.datatype=16;
                    %niiheader.hdr.dime.bitpix=32;
                    outfilename = [T1MAPS_folder,slsh,the_filename,'_mFAfit_T1MAPS'];
                    %outfilename = [T1MAPS_folder,slsh,'T1MAPS',slsh,the_filename,'_mFAfit_T1MAPS'];
                    sdir=fileparts(outfilename);
                    mkdir(sdir);
                    
                    
                    if use_b1_in_fit==1
                        if DO_SMOOTH_B1
                            if EXTREME_SMOOTH_B1
                                outfilename=char([outfilename,'_EXTREMEsmoothB1']);
                            else
                                outfilename=char([outfilename,'_smoothB1']);
                            end
                        else
                            outfilename=char([outfilename,'_withB1']);
                        end
                    end
            end
            
            
            %outfilename=[outfilename '_' mydate];
            %save_untouch_nii(mynii,outfilename);
            nii_tool('save',mynii,[outfilename,'.nii.gz']);
            
%             DO_SMOOTH_B1=1
%             EXTREME_SMOOTH_B1=0
%             DO_NOT_USE_B1_DATA=0
             

        end
        
        
        %% Save matlab data to show fits
        if ~exist('dataB1','var')
            dataB1=0;
        end
        
        if ~exist('FAmaps','var')
            FAmaps=0;
        end
        
        if ~exist('TR','var')
            TR=0;
        end
        
        if ~exist('TIs','var')
            TIs=0;
        end
        save([outfilename,'.mat'],'T1maps','FAmaps','M0maps','TIs','dataT1','dataB1', 'use_b1_in_fit','fit_type','R2adjusted','fitparam','TR')
        %saveinparfor([outfilename,'.mat'], T1maps,FAmaps,M0maps,TIs,dataT1,dataB1, use_b1_in_fit,fit_type,R2adjusted',fitparam,TR);
        
        %% TicTocs
         disp(datetime);
        patient_time=num2str(round(toc(patienttic)/60,1));
        disp(['Time required to process ', pfolders{patient} ,' was ',patient_time ,' minutes',]);
        
        diary off
        processed=processed+1;
        
        
end



 disp([' ']);
    disp(['Finished! I managed to process ', num2str(processed),' out of the ' ,num2str(numel(pfolders)),' datasets selected. ']);
    disp(pfolders(:));
    batchtime=num2str(round(toc(batchtic)/60,1));
    
    disp(['It took ',batchtime ,' minutes to run them.',]);