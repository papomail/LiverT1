%% T1 maps calculation

%set FSL invironment
setenv('FSLDIR','/usr/local/fsl/')
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
fclose all
format long
clear all

pathtoadd=genpath(pwd);
addpath(pathtoadd)
%% Options
SAVE_NII=1
MOCO=1

%% Input parameters

TIs_CEST = [100 235 550 1300 3000 7000]; %ms
TIs_MOLLI3 = [120, 200, 280]; %ms
TIs_MOLLI5 = [120, 200, 280, 360, 420]; %ms
TIs_LIVER = [100 200 300 500 1000 2000]; %ms

TIs = TIs_LIVER;

%FA = 7; 
%TR = 6000; %ms


FAs=[2.5 5 7 10 15];

%% Set original path
[the_folder] = pathsetting('Desktop'); 

%%  MOCO
MOCO=0
if MOCO==1
    tic
    [myfileselection]=MOCO_LIVER_T1;%(myfileselection);
    toc
end
    
%% Select the data
message = 'Please select folder with the dicoms OR the inversion-recovery niftis ';
[data_cell, folderfiles, folders] = reader_Aug2018(the_folder,message);
num_dataSets = numel(data_cell);
%%Read the selected data
clear dataT1 dataT1_init sequenceTI sequenceFA


%%
% Load and check if InversionRecovery of MultiFA fit
NIFTIfolder=fileparts(folderfiles{1}{1});
dicomheaders=char([NIFTIfolder,'/dcmHeaders.mat']);

load(dicomheaders)


for ii=1:num_dataSets
    
    fullpath=[data_cell{ii}.fileprefix,'.nii'];
    [~,myfilename,~]=fileparts(fullpath);
    
    myfs=strfind(myfilename,'_FLIRT');  %% Find the headers of the original (non-registered) file.
    if myfs
        myfilename=myfilename(1:myfs-1);
            [~,myfilename,~]=fileparts(myfilename);
    end
    
    if isfield(h.(myfilename), 'InversionTime')
        sequenceTI(ii)=h.(myfilename).InversionTime;
        fit_type= 'inversion_recovery';
        LookLocker='no'
        mD_IR='no'
    elseif isfield(h.(myfilename), 'TriggerTime'); %% For MOLLI/LookLocker
        sequenceTIint=220 ;
        sequenceTI1=h.(myfilename).TriggerTime;
        sequenceTInum=h.(myfilename).NumberOfPhasesMR;
        sequenceTI=double([sequenceTI1:sequenceTIint:sequenceTInum*sequenceTIint])
        disp('Please check if these TIs are correct')
        fit_type= 'inversion_recovery';
        LookLocker='yes'
        mD_IR='no'
    elseif strncmp(h.(myfilename).NiftiName, 'mD_RPP_TI',9); %% For MOLLI/LookLocker
        
        myTI1=str2double(h.(myfilename).NiftiName(10:13));
        if isnan(myTI1)
            myTI1=str2double(h.(myfilename).NiftiName(10:12));
        end
        sequenceTI(ii)=myTI1;
        fit_type= 'inversion_recovery';
        mD_IR='yes'
        LookLocker='no'
    else
        sequenceFA(ii)=h.(myfilename).FlipAngle;
        TR(ii)=h.(myfilename).RepetitionTime;
        fit_type= 'multi_fa';
    end
    
    
    
    %%%% Scaling .nii for matlab to display float point
    niiSS = double(data_cell{ii}.hdr.dime.scl_slope);
    %dicomSS= double(h.(myfilename).MRScaleSlope);
    %dicomRS= double(h.(myfilename).RescaleSlope);
        

    niiRI= double(data_cell{ii}.hdr.dime.scl_inter);
   % dicomSI= double(h.(myfilename).MRScaleIntercept);
   % dicomRI= double(h.(myfilename).RescaleIntercept);

    
    
%     if niiSS==0 || dicomSS==0%%(as in Siemens)
%         niiSS=1;
%         dicomSS=1;
%         niiRI=0;
%         dicomRI=0;
%         disp('No Scaling Slope was found. Using SS=1')
%     end
    
    if 0%myfs
     dataT1_init(:,:,:,:,ii)= double(data_cell{ii}.img)/dicomRS+dicomRI/(dicomRI*dicomRS) ;
    elseif 1
    dataT1_init(:,:,:,:,ii)= double(data_cell{ii}.img)*niiSS+niiRI;
    elseif 0
    dataT1_init(:,:,:,:,ii)= double(data_cell{ii}.img)/niiSS/dicomSS+dicomRI/(dicomRS*dicomSS) ;
    elseif 0
            dataT1_init(:,:,:,:,ii)= double(data_cell{ii}.img)/niiSS/dicomSS*dicomRS; % BE CAREFULL!!!!
    end
    
    %dicomRS_A(ii)=dicomRS;
    %dicomSS_A(ii)=dicomSS;
    %niiSS_A(ii)=niiSS;
   
    % Scale to Floating Point: (Philips) (and Torbens's slides)
    % FP= SV/SS + RI/(RS*SS), where FP: "Floating Point"  
  
    
    %QuantitativeValue = [SV-ScInt] / ScSlp,
%   where SV = stored DICOM pixel value, ScInt = Scale Intercept = (2005,100d), ScSlp = Scale Slope = (2005,100e).
    
    
    % SV: "Stored Value"; stored in DICOM image file
    % SS: "MR Scale Slope"; DICOM header tag: (2005, 100E)
    % RI: "Rescale Intercept"; DICOM header tag: (0028, 1052)
    % RS: "Rescale Slope"; DICOM header tag: (0028, 1053)
    % These values have been scaled already.
end

if strcmp(fit_type , 'multi_fa')
    if mean(TR)==TR(1) %%chek if all TR are the same
        TR=double(TR(1));
    else
        disp({'WHAT! Your TR values are not the same!'; 'What are you trying to accomplish with this?? ';'It wont work...'});
        
    end
end



%% Sort and Reorient
switch fit_type
    case 'inversion_recovery'
        
        [newTIs,perTIs]=sort(sequenceTI);
        TIs=newTIs;
        
        if strcmp(LookLocker,'yes' )
            reorientdataT1(:,:,1,1,:)=dataT1_init;%If it's LLocker, is one slice only and alredy sorted by TI (TriggerTime)
        else
            
            reorientdataT1=dataT1_init(:,:,:,:,perTIs);%sort based on TI
        end
        % reorient image
        reorientdataT1=flip(rot90(reorientdataT1),2);%% rotate and flip image
        
        if strcmp(mD_IR , 'yes')
            showFatWater(reorientdataT1)% show 4 images (fat water outofphase inphase)

reorientdataT1=flip(reorientdataT1(:,:,:,[1,2],:),4);  %% Put WaterOnly first and FatOnly second, and get rid of the rest.
%reorientdataT1=flip(reorientdataT1(:,:,:,[3,4],:),4);  %% Put INPHASE first and OUTOFPHASE second, and get rid of the rest.

showFatWater(reorientdataT1)% show 2 images (water fat) 
        end
        
    case 'multi_fa'
        [newFAs,perFAs]=sort(sequenceFA);
        FAs=newFAs;
        if FAs(1)==2
            FAs(1)=2.5;  %%% BECAUSE THE DICOM HEADER SAYS 2 DEGREES INSTEAD OF 2.5!!
        end
        reorientdataT1=dataT1_init(:,:,:,:,perFAs); %sort based on FA
% reorient image
reorientdataT1=permute(reorientdataT1,[1 3 2 4 5]);%% permute dimensions 2 and 3 to put images as yxz
reorientdataT1=rot90(reorientdataT1); %% rotate image
showFatWater(reorientdataT1)% show 4 images (fat water outofphase inphase)  (ORDER IN IFALD2: f i o w )

reorientdataT1=flip(reorientdataT1(:,:,:,[1,2],:),4);  %% Put WaterOnly first and FatOnly second, and get rid of the rest.
%reorientdataT1=flip(reorientdataT1(:,:,:,[3,4],:),4);  %% Put INPHASE first and OUTOFPHASE second, and get rid of the rest.

showFatWater(reorientdataT1)% show 2 images (water fat) 

end


%% Load B1map if aquired.
B1map_file=char([NIFTIfolder,'/B1map_Liver_coronal.nii']);
if ~(exist(B1map_file,'file')==2)
B1map_file=char([NIFTIfolder,'/WIP_B1map_Liver_coronal.nii']);
 ICanFindB0Map=(exist(B1map_file,'file')==2);
else
   ICanFindB0Map=1;
end

if ICanFindB0Map
    disp('Loading B1map.')
    B1map_data=load_untouch_nii(B1map_file);
    B1map_data=B1map_data.img;
    
    B1map_data=permute(B1map_data,[1 3 2 4 5]);%% reorient B1maps
    B1map_data=rot90(B1map_data);
    B1map_data=B1map_data(:,:,:,1);
    use_b1_in_fit=1;
    % interpolate B1map to match mFA data
    disp('Will assume the same FOV as the data for T1 mapping.')
    
    dvsize=size(reorientdataT1);
    xq=linspace(1,size(B1map_data,1),dvsize(1));
    yq=linspace(1,size(B1map_data,2),dvsize(2));
    zq=linspace(1,size(B1map_data,3),dvsize(3));
    [Xq,Yq,Zq]=ndgrid(xq,yq,zq);
    
    xo=1:size(B1map_data,1);
    yo=1:size(B1map_data,2);
    zo=1:size(B1map_data,3);
    [Xo,Yo,Zo]=ndgrid(xo,yo,zo);
    
    B1map_data_int=interpn(Xo,Yo,Zo,B1map_data,Xq,Yq,Zq);
else
        use_b1_in_fit=0;

    warning('I could not find any B1 maps. Continuing without it.')
    B1map_data_int=ones(size(reorientdataT1,1),size(reorientdataT1,2),size(reorientdataT1,3));
end


%% select ROI in image

figure, imagesc(reorientdataT1(:,:,round(size(reorientdataT1,3)/2),1,1))
title('Please draw a rectangle in the liver where fitting is required. [double clic to finish]')
h = imrect;
position = round(wait(h));
close
xrange=position(1):position(1)+position(3);
yrange=position(2):position(2)+position(4);
reorientdataT1_s=reorientdataT1(yrange,xrange,:,:,:);

%figure, imagesc(reorientdataT1_s(:,:,round(size(reorientdataT1_s,3)/2),1,1))

slrange=':';

switch fit_type
    case 'multi_fa'
        slrange=[26];%[51:1:75] % Take the useful middle slices
        dataT1=reorientdataT1_s(:,:,slrange,1,:); %% fit the water signal only
        dataB1=B1map_data_int(yrange,xrange,slrange,:,:);
        
        mergezslizes=0%% merges zslices in pairs to obtain similar slice thickness as in the IR_TSE sequence
        if mergezslizes
        dataT1odd=dataT1(:,:,1:2:end,:,:);%% merge zslices in pairs
        dataT1even=dataT1(:,:,2:2:end,:,:);
        dataT1merg=(dataT1odd+dataT1even)/2;
        dataT1=dataT1merg;
        
         dataB1odd=dataB1(:,:,1:2:end);%% merge zslices in pairs
        dataB1even=dataB1(:,:,2:2:end);
        dataB1merg=(dataB1odd+dataB1even)/2;
        dataB1=dataB1merg;
        clear  dataT1_mer dataB1merg
        end
        
    case 'inversion_recovery'
        if strcmp(LookLocker,'yes')
            clear dataT1
        dataT1(:,:,1,1,:)=reorientdataT1_s;
        slrange=[1]
                dataB1=B1map_data_int(yrange,xrange,slrange,:,:);

        else
        slrange=[4:17] % Take the useful middle slices
                slrange=[1] % Take the useful middle slices
                
dataT1=reorientdataT1_s(:,:,slrange,1,:); %% fit the water signal only
                dataB1=B1map_data_int(yrange,xrange,slrange,:,:);

        end
end

 imagesize=size(dataT1);
 imagesize=imagesize(1:end-1);
% new path
filepath_T1data=cell2mat(folderfiles{1}(1)); 

%  check voxel(s) plot(s)
fatorwat=1;% in multiFA data: 1 for water; 2 for fat
sl=round(size(dataT1,3)/2);
figure,imagesc(dataT1(:,:,sl,fatorwat,1))
title('Choose the few voxels to check the raw data. [double clic to finish]')
[pxs, pys] = getpts;
close
pxs=round(pxs);
pys=round(pys);
for ii=1:numel(pxs)
figure,plot(squeeze(dataT1(pys(ii),pxs(ii),sl,fatorwat,:)),'-o')
title('Data in the selected voxel')
end


%% Fitting
switch fit_type
    case 'inversion_recovery'
        disp('starting IR fit')
        tic
        if use_b1_in_fit
           [fitparam]= IR_T1FITTING_Oct2018(dataT1,TIs(:),dataB1);
        %[fitparam]=T1FITTING_Pat_Aug2018(dataT1,TIs(:),dataB1);
        else
           [fitparam]= IR_T1FITTING_Oct2018(dataT1,TIs(:));
        %[fitparam]=T1FITTING_Pat_Aug2018(dataT1,TIs(:));
        end
        toc
        xvals=TIs;
    case 'multi_fa'
        disp('starting Multi FA fit')
        tic
        if use_b1_in_fit
            [fitparam] = multiFAfit(dataT1,FAs,TR,dataB1); %fit with B1
            
        else
            [fitparam] = multiFAfit(dataT1,FAs,TR);%fit without B1
        end
        toc
        xvals=FAs;
end
  

 %% R^2 calculation
myfval=fitparam.fval;
[ R2adjusted ] = R2calc( dataT1,xvals,myfval );
%% mask out data that fit not a good fit.
do_mask_based_on_R2=1
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



%% 
switch fit_type
    case 'inversion_recovery'
        if use_b1_in_fit==0
            FlipAngles=fitparam.FlipAngle.*rsquarefilt;
            FAmaps=reshape(FlipAngles,size(T1s,1),size(T1s,2)*size(T1s,3));
        end
        
        if strcmp(LookLocker,'yes')
            myA=M0s;
            myB=M0s.*(1+FlipAngles);
            LLockerFactor=(myB./myA-1);
            T1LLfit= T1s;
            T1s=T1s.*LLockerFactor;
        end
        
end

%% reshape for plots 
T1maps=reshape(T1s,size(T1s,1),size(T1s,2)*size(T1s,3));
M0maps=reshape(M0s,size(T1s,1),size(T1s,2)*size(T1s,3));

    
%aamaps=reshape(fitparam.aa,size(fitparam.RelaxTime,1),size(fitparam.RelaxTime,2)*size(fitparam.RelaxTime,3),size(fitparam.RelaxTime,4));


% plot T1maps
caxmax=2;
caxmin=0;



figure,
imagesc((T1maps/1000)),title({'T1 map [s]',['(Masked with Adjusted R^2 > ',num2str(myr2lim),')']}), colorbar
colormap(jet)
caxis([caxmin caxmax])

% 
% %plot M0maps
% figure,
% imagesc(M0maps),title({'M0 map [s]',['(Masked with Adjusted R^2 > ',num2str(myr2lim),')']}), colorbar
% colormap(jet)
% caxis([-0 400])

if use_b1_in_fit==0
%plot FA
myFlipA=180/pi*real(acos(-FAmaps));
figure,
imagesc(myFlipA),title({'Flip angle [degrees]',['(Masked with Adjusted R^2 > ',num2str(myr2lim),')']}), colorbar
colormap(jet)
caxis([90 180])
end

  %% SAVE .nii MAPS
if SAVE_NII==1
[the_folder, the_filename] = fileparts(char(filepath_T1data));


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
 mynii=data_cell{1};
 
 %get the rid of scaling slope
 mynii.hdr.dime.scl_slope=1;          
 mynii.hdr.dime.scl_inter=0;
 
 switch fit_type
     case 'inversion_recovery'
         fullsizeim=flip(rot90(double(data_cell{1}.img)),2); %open original image to substitute with t1maps values (xray effect)
         myint=fullsizeim(:,:,round(size(fullsizeim,3)/2),1,1);
         myint(myint==0)=nan;
         fullsizeim=fullsizeim*200/nanmedian(myint(:)); %fill the original image with t1maps values (xray effect)
         fullsizeim(yrange,xrange,slrange)=tosave_T1maps;
         fullsizeim=rot90(flip(fullsizeim,2),-1); %Re-orient to initial orientation
         mynii.img=fullsizeim;
        
         
         
         if strcmp(LookLocker,'yes')
             outfilename = [the_folder,slsh,the_filename,'_LLfit_T1MAPS'];
             
         else
             outfilename = [the_folder,slsh,the_filename,'_IRfit_T1MAPS'];
         end
         
          if use_b1_in_fit==1
         outfilename=char([outfilename,'_withB1']);
          end
         
     case 'multi_fa'
         
         fullsizeim=rot90(permute(double(data_cell{1}.img),[1 3 2 4])); %open original image to substitute with t1maps values (xray effect)
         fullsizeim=fullsizeim(:,:,:,4);%take the water only image
         myint=fullsizeim(:,:,round(size(fullsizeim,3)/2),1,1);
         myint(myint==0)=nan;
         fullsizeim=fullsizeim*200/nanmedian(myint(:)); %fill the original image with t1maps values (xray effect)
         
         
         if mergezslizes
             tosave_T1maps2=zeros(size(tosave_T1maps,1),size(tosave_T1maps,2),size(tosave_T1maps,3)*2); %put every fittes slice into 2 original-thinkness slices
             tosave_T1maps2(:,:,1:2:end)=tosave_T1maps;
             tosave_T1maps2(:,:,2:2:end)=tosave_T1maps;
             
             fullsizeim(yrange,xrange,slrange)=tosave_T1maps2;
         else
             fullsizeim(yrange,xrange,slrange)=tosave_T1maps;
         end
         
         fullsizeim=permute(rot90(fullsizeim,-1),[1 3 2 4]); %Re-orient to initial orientation
         mynii.img=fullsizeim;
         
         
         
         %make sure to go one dimension less as we now have maps
         mynii.hdr.dime.dim(2:5) = [size(fullsizeim),1];
         %mynii.hdr.dime.pixdim(5) = 0;
         
         %make sure output is float
         %niiheader.hdr.dime.datatype=16;
         %niiheader.hdr.dime.bitpix=32;
         outfilename = [the_folder,slsh,the_filename,'_mFAfit_T1MAPS'];
         if use_b1_in_fit==1
         outfilename=char([outfilename,'_withB1']);
         end
 end
 
 save_untouch_nii(mynii,outfilename);
 
 
 

end 

%% See the fits 
if exist('FAmaps','var')
    PlotT1Fits(T1maps,FAmaps,M0maps,TIs,dataT1,1, use_b1_in_fit,fit_type,R2adjusted)

else
    PlotT1Fits(T1maps,1,M0maps,TIs,dataT1,dataB1, use_b1_in_fit,fit_type,R2adjusted)
end

%% Select voxel to plot
figure,

imagesc(T1maps/1000),title('Select the voxels you want to see the fit on. (Double-click to finish)'), colorbar
caxis([caxmin caxmax])
colormap(jet)
set(gca,'XTick','')
set(gca,'YTick','')
%xticks('')
%yticks('')
if isfield(fitparam,'Aconstant')
                Aconstant=fitparam.Aconstant;
                Aconstantmap=reshape(Aconstant,size(Aconstant,1),size(Aconstant,2)*size(Aconstant,3),size(Aconstant,4));
end
[pxs, pys] = getpts;

%Plot selected voxels  (only works for IR for the moment)
for ii=1:numel(pxs)
    px=round(pxs(ii));
    py=round(pys(ii));
    
    %aa1=squeeze(aamaps(py,px,:,:));
    M01=squeeze(M0maps(py,px,:,:));
    if exist('LookLocker','var') && strcmp(LookLocker,'yes')
        T11=squeeze(T1LLfit(py,px,:,:));
    else
        T11=squeeze(T1maps(py,px,:,:));
    end
    
    reshaped_dataT1=reshape(dataT1,size(dataT1,1),size(dataT1,2)*size(dataT1,3),size(dataT1,4),size(dataT1,5));
    dataT11=squeeze(reshaped_dataT1(py,px,:,:));
    
    
    switch fit_type
        case 'inversion_recovery'
            
            figure,
            plot(TIs,squeeze(dataT11(:,1)),'o')
            hold on
            %plot(myxx,aa1(1)+abs(M01(1)*exp(-myxx/T11(1))+Minf1(1)*(1-exp(-myxx/T11(1) ))))
            
            if use_b1_in_fit
                B1map=reshape(dataB1,size(dataB1,1),size(dataB1,2)*size(dataB1,3));
                Minf1=-cos(pi*B1map(py,px,:,:)/100);
            else
                Minf1=squeeze(FAmaps(py,px,:,:));
            end
            
            myxx=logspace(0,3.8,300);
            plot(myxx,M01*abs(1-(1+Minf1)*exp(-myxx/T11)))
            xlabel('Inversion delay (ms)')
            ylabel('Signal Intensity')
        case 'multi_fa'
            
            reshaped_FAcor=reshape(fitparam.FA_b1corrected,size(dataT1,1),size(dataT1,2)*size(dataT1,3),size(dataT1,4),size(dataT1,5));
            FAcor1=squeeze(reshaped_FAcor(py,px,:,:));
            
            myxx=linspace(0,20,300);
            figure,
            plot(FAcor1,squeeze(dataT11(:,1)),'o')
            hold on
            
            if exist('Aconstantmap','var')
                Aconstantmap1=squeeze(Aconstantmap(py,px,:,:));
            else
                Aconstantmap1=0;
            end
            plot(myxx,Aconstantmap1+M01*sin(myxx*pi/180).*(1-exp(-TR/T11))./(1-cos(myxx*pi/180).*exp(-TR/T11)))
            xlabel('Efective FlipAngle [degrees]')
            ylabel('Signal Intensity')
    end
    
end



STOP

%% Registration bit (incomplete)
DO_REGISTRATION=0;
if DO_REGISTRATION
    registration_Aug2018
end




%% Select ROI for single fit
figure,

    imagesc(T1maps/1000),title('Select ROI to do an extra fit on. (Double-click to finish)'), colorbar
    caxis([caxmin caxmax])
colormap(jet)
   set(gca,'XTick','')
   set(gca,'YTick','')
%xticks('')
%yticks('')

h = imrect;
position = round(wait(h));
close
xrange=position(1):position(1)+position(3);
yrange=position(2):position(2)+position(4);



rsquarefilt_mat=reshape(rsquarefilt,[size(rsquarefilt,1),size(rsquarefilt,2)*size(rsquarefilt,3)]); 
rsquarefilt_mat=repmat(rsquarefilt_mat,[1,1,size(reshaped_dataT1,4)]);

%dataT1_ROI=squeeze(reshaped_dataT1(yrange,xrange,:,:)).*rsquarefilt_mat(yrange,xrange,:,:);%% Mask out voxels with R^2<09 fot the single fit
dataT1_ROI=squeeze(reshaped_dataT1(yrange,xrange,:,:));%% 

[T1cal ,T1conf_low ,T1conf_high]=ROI_Fit(TIs,dataT1_ROI);

%% plots
clow=0;
chigh=2000;

    figure,
    colormap(jet)
   
    subplot(2,3,2)
    imagesc(T1cal)
    xlabel('T1 [ms]')
    caxis([clow chigh])
         set(gca,'XTick','')
         set(gca,'YTick','')
         colorbar
    
    subplot(2,3,1)
    imagesc(T1conf_low)
    xlabel( 'Lower bound T1 (with 95% confidence) [ms]')
    caxis([clow chigh])
         set(gca,'XTick','')
         set(gca,'YTick','')
         colorbar
    
    subplot(2,3,3)
    imagesc(T1conf_high),
    xlabel('Upper bound T1 (with 95% confidence) [ms]')
    caxis([clow chigh])
         set(gca,'XTick','')
         set(gca,'YTick','')
         colorbar
         
          subplot(2,3,6)
    histogram(T1conf_high),
    xlabel('Histogram Upper bound T1')
    caxis([clow chigh])
    
     subplot(2,3,5)
    histogram(T1cal),
    xlabel('Histogram T1')
    caxis([clow chigh])
        
     subplot(2,3,4)
    histogram(T1conf_low),
    xlabel('Histogram Lower bound T1')
    caxis([clow chigh])
  
filename = [the_folder,slsh,the_filename,'_HISTO_T1MAPS'];
saveas(gcf,filename)      
         

      figure,imagesc(T1conf_high - T1conf_low),
    title('\DeltaT1 within 95% confidence [ms]' )
    colormap(jet)
    caxis([0 2000])
      set(gca,'XTick','')
         set(gca,'YTick','')
         colorbar
    
    
nanmedian(T1cal(:))
nanmean(T1cal(:))

filename = [the_folder,slsh,the_filename,'_DELTA_T1MAPS'];
saveas(gcf,filename)

%% remove newly added path
rmpath(pathtoadd)