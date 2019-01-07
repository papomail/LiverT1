% T1 maps calculation
%% 
fclose all
format long
clear all

pathtoadd=genpath(pwd);
addpath(pathtoadd)
%% Options
SAVE_NII=1


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
message = 'Please select folder with the dicoms OR the inversion-recovery niftis ';

%% Select the data
[data_cell, folderfiles, folders] = reader_Aug2018(the_folder,message);
num_dataSets = numel(data_cell);
%%Read the selected data
clear dataT1 dataT1_init sequenceTI sequenceFA 



% Load and check if InversionRecovery of MultiFA fit
NIFTIfolder=fileparts(folderfiles{1}{1});
dicomheaders=char([NIFTIfolder,'/dcmHeaders.mat']);

load(dicomheaders)
 

for ii=1:num_dataSets 
   
     fullpath=[data_cell{ii}.fileprefix,'.nii'];
    [~,myfilename,~]=fileparts(fullpath);
    if isfield(h.(myfilename), 'InversionTime')
        sequenceTI(ii)=h.(myfilename).InversionTime;
        fit_type= 'inversion_recovery';
        LookLocker='no'
    elseif isfield(h.(myfilename), 'TriggerTime'); %% For MOLLI/LookLocker
         sequenceTIint=220 ;
         sequenceTI1=h.(myfilename).TriggerTime;
         sequenceTInum=h.(myfilename).NumberOfPhasesMR;
         sequenceTI=double([sequenceTI1:sequenceTIint:sequenceTInum*sequenceTIint])
         disp('Please check if these TIs are correct')
        fit_type= 'inversion_recovery';
        LookLocker='yes'
    else
        sequenceFA(ii)=h.(myfilename).FlipAngle;
        TR(ii)=h.(myfilename).RepetitionTime;
        fit_type= 'multi_fa';
    end
    
        
    
    %%%% Scaling .nii for matlab to display float point
     %SS = data_cell{ii}.hdr.dime.scl_slope;   
      SS= double(h.(myfilename).MRScaleSlope);
     %RI= data_cell{ii}.hdr.dime.scl_inter;
        RI= double(h.(myfilename).MRScaleIntercept);

     if SS==0 %%(Siemens)
         SS=1;
         RI=0;
         disp('No Scaling Slope was found. Using SS=1')
     end
    
             dataT1_init(:,:,:,:,ii)= double(data_cell{ii}.img)/SS+RI/(RI+SS) ;
             
% Scale to Floating Point: (Philips)
% FP= SV/SS + RI/(RS*SS), where FP: "Floating Point"

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
        
    case 'multi_fa'
        [newFAs,perFAs]=sort(sequenceFA);
        FAs=newFAs;
        reorientdataT1=dataT1_init(:,:,:,:,perFAs); %sort based on FA
% reorient image
reorientdataT1=permute(reorientdataT1,[1 3 2 4 5]);%% permute dimensions 2 and 3 to put images as yxz
reorientdataT1=rot90(reorientdataT1); %% rotate image
showFatWater(reorientdataT1)% show 4 images (fat water outofphase inphase)

reorientdataT1=flip(reorientdataT1(:,:,:,[1,4],:),4);  %% Put WaterOnly first and FatOnly second, and get rid of the rest.
showFatWater(reorientdataT1)% show 2 images (water fat) 


        % Load B1map if aquired/needed.
        B1map_file=char([NIFTIfolder,'/WIP_B1map_Liver_coronal.nii']);
        B1map_data=load_untouch_nii(B1map_file);

        if isempty(B1map_data)
            warning('I could not find any B1 maps. Continuing without it.')
        else
            disp('Loading B1map.')
            B1map_data=B1map_data.img;
        end
        B1map_data=permute(B1map_data,[1 3 2 4 5]);%% reorient B1maps
        B1map_data=rot90(B1map_data);
        B1map_data=B1map_data(:,:,:,1);
        
        % interpolate B1map to match mFA data
       disp('Will assume the same FOV as the multiFA data.')

       dvsize=size(reorientdataT1);
       xq=linspace(1,size(B1map_data,1),dvsize(1));
       yq=linspace(1,size(B1map_data,2),dvsize(2));
       zq=linspace(1,size(B1map_data,3),dvsize(3));
       tgrid=meshgrid(xq,yq,zq);
       
       xo=1:size(B1map_data,1);
       yo=1:size(B1map_data,2);
       zo=1:size(B1map_data,3);
       ogrid=meshgrid(xo,yo,zo);
       
B1map_data_int=interpn(ogrid,B1map_data,tgrid);

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
        slrange=[41:70] % Take the useful middle slices
        dataT1=reorientdataT1_s(:,:,slrange,1,:); %% fit the water signal only
        dataT1odd=dataT1(:,:,1:2:end,:,:);%% merge zslices in pairs
        dataT1even=dataT1(:,:,2:2:end,:,:);
        dataT1merg=(dataT1odd+dataT1even)/2;
        dataT1=dataT1merg;
        clear  dataT1_mer
    case 'inversion_recovery'
        if strcmp(LookLocker,'yes')
            clear dataT1
        dataT1(:,:,1,1,:)=reorientdataT1_s;
        slrange=[1]
        else
        slrange=[10:10] % Take the useful middle slices
        dataT1=reorientdataT1_s(:,:,slrange,1,:); %% fit the water signal only
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
        [fitparam]=T1FITTING_Pat_Aug2018(dataT1,TIs(:));
        toc;
        xvals=TIs;
    case 'multi_fa'
        disp('starting Multi FA fit')
        tic
        [fitparam] = multiFAfit(dataT1,FAs,TR);
        toc 
        xvals=FAs;
end
  

 %% R^2 calculation
myfval=fitparam.fval;
[ R2adjusted ] = R2calc( dataT1,xvals,myfval );
%% mask out data that fit not a good fit.
myr2lim=0.0;
rsquarefilt=R2adjusted>myr2lim;
rsquarefilt=double(rsquarefilt);
rsquarefilt(rsquarefilt==0)=nan;


T1s=fitparam.RelaxTime.*rsquarefilt;
M0s=fitparam.M0.*rsquarefilt;
%myfval=myfval.*rsquarefilt;



%% 
switch fit_type
    case 'inversion_recovery'
FlipAngles=fitparam.FlipAngle.*rsquarefilt;
FAmaps=reshape(FlipAngles,size(T1s,1),size(T1s,2)*size(T1s,3));


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
caxmax=3;
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
% 
% %plot FA
% myFlipA=180/pi*real(acos(-FAmaps));
% figure,
% imagesc(myFlipA),title({'Flip angle [degrees]',['(Masked with Adjusted R^2 > ',num2str(myr2lim),')']}), colorbar
% colormap(jet)
% caxis([90 180])


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
  fullsizeim=fullsizeim*0.5; %fill the original image with t1maps values (xray effect)
  fullsizeim(yrange,xrange,slrange)=tosave_T1maps;
  fullsizeim=rot90(flip(fullsizeim,2),-1); %Re-orient to initial orientation
  mynii.img=fullsizeim;
  if strcmp(LookLocker,'yes')
        outfilename = [the_folder,slsh,the_filename,'_LLfit_T1MAPS'];

  else
  outfilename = [the_folder,slsh,the_filename,'_IRfit_T1MAPS'];
  end
     case 'multi_fa'
      
 fullsizeim=rot90(permute(double(data_cell{1}.img),[1 3 2 4])); %open original image to substitute with t1maps values (xray effect)
 fullsizeim=fullsizeim(:,:,:,4);%take the water only image
 fullsizeim=fullsizeim*0.5; %fill the original image with t1maps values (xray effect)
   
   tosave_T1maps2=zeros(size(tosave_T1maps,1),size(tosave_T1maps,2),size(tosave_T1maps,3)*2); %put every fittes slice into 2 original-thinkness slices
  tosave_T1maps2(:,:,1:2:end)=tosave_T1maps;
    tosave_T1maps2(:,:,2:2:end)=tosave_T1maps;

   fullsizeim(yrange,xrange,slrange)=tosave_T1maps2;
       fullsizeim=permute(rot90(fullsizeim,-1),[1 3 2 4]); %Re-orient to initial orientation
 mynii.img=fullsizeim;

 
 
 %make sure to go one dimension less as we now have maps
 mynii.hdr.dime.dim(2:5) = [size(fullsizeim),1];
 %mynii.hdr.dime.pixdim(5) = 0;

 %make sure output is float
 %niiheader.hdr.dime.datatype=16;
 %niiheader.hdr.dime.bitpix=32;
          outfilename = [the_folder,slsh,the_filename,'_mFAfit_T1MAPS'];

 end
 
 save_untouch_nii(mynii,outfilename);
 
 
 

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

[pxs, pys] = getpts;


%Plot selected voxels
for ii=1:numel(pxs)
    px=round(pxs(ii));
py=round(pys(ii));

%aa1=squeeze(aamaps(py,px,:,:));
M01=squeeze(M0maps(py,px,:,:));
 if strcmp(LookLocker,'yes')
     T11=squeeze(T1LLfit(py,px,:,:));
 else
T11=squeeze(T1maps(py,px,:,:));
 end
Minf1=squeeze(FAmaps(py,px,:,:));

reshaped_dataT1=reshape(dataT1,size(dataT1,1),size(dataT1,2)*size(dataT1,3),size(dataT1,4),size(dataT1,5));
dataT11=squeeze(reshaped_dataT1(py,px,:,:));
myxx=logspace(0,3.8,300);


    figure,
    plot(TIs,squeeze(dataT11(:,1)),'o')
    hold on
    %plot(myxx,aa1(1)+abs(M01(1)*exp(-myxx/T11(1))+Minf1(1)*(1-exp(-myxx/T11(1) ))))

    
plot(myxx,M01(1)*abs(1-(1+Minf1(1))*exp(-myxx/T11(1))))


end




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

dataT1_ROI=squeeze(reshaped_dataT1(yrange,xrange,:,:)).*rsquarefilt_mat(yrange,xrange,:,:);%% Mask out voxels with R^2<09 fot the single fit

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
        
         

      figure,imagesc(T1conf_high - T1conf_low),
    title('\DeltaT1 within 95% confidence [ms]' )
    colormap(jet)
    caxis([0 2000])
      set(gca,'XTick','')
         set(gca,'YTick','')
         colorbar
    
    
nanmedian(T1cal(:))
nanmean(T1cal(:))




%% remove newly added path
rmpath(pathtoadd)