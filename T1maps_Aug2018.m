% T1 maps calculation
%% 
fclose all
format long
clear all


addpath(genpath(pwd))
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
%% Set original path
[the_folder] = pathsetting('Desktop'); 
message = 'Please select folder with the dicoms OR the inversion-recovery niftis ';

%% Select the data
[data_cell, folderfiles, folders] = reader_Aug2018(the_folder,message);
num_dataSets = numel(data_cell);
%%Read the selected data
clear dataT1 dataT1_init
if num_dataSets == 2
    data_cell{1} = imresize(data_cell{2},[size(data_cell{2},1),size(data_cell{2},2)]);
end
for ii=1:num_dataSets 
   
    
    %%%% Scaling .nii for matlab to display float point
     SS = data_cell{ii}.hdr.dime.scl_slope;          
     RI= data_cell{ii}.hdr.dime.scl_inter;

     if SS==0 %%(Siemens)
         SS=1;
         RI=0;
         disp('No Scaling Slope was found. Using SS=1')
     end
    
             dataT1_init(:,:,:,:,ii)= double(data_cell{ii}.img)/SS+RI/(RI+SS);
             
% Scale to Floating Point: (Philips)
% FP= SV/SS + RI/(RS*SS), where FP: "Floating Point"

% SV: "Stored Value"; stored in DICOM image file
% SS: "MR Scale Slope"; DICOM header tag: (2005, 100E)
% RI: "Rescale Intercept"; DICOM header tag: (0028, 1052)
% RS: "Rescale Slope"; DICOM header tag: (0028, 1053)
 % These values have been scaled already. 
end


%% SORT images based on their TI time
NIFTIfolder=fileparts(folders{1});
dicomheaders=char([NIFTIfolder,'/dcmHeaders.mat']);
load(dicomheaders)
for ii=1:num_dataSets
    fullpath=[data_cell{ii}.fileprefix,'.nii'];
    [~,myfilename,~]=fileparts(fullpath);
    sequenceTI(ii)=h.(myfilename).InversionTime;
end

[newTIs,perTIs]=sort(sequenceTI);
dataT1_init(:,:,:,:,1:num_dataSets)=dataT1_init(:,:,:,:,perTIs);
TIs=newTIs; 

%% reorient image

reorientdataT1=flip(rot90(dataT1_init),2);
%% select ROI in image

figure, imagesc(reorientdataT1(:,:,round(size(reorientdataT1,3)/2),1,1))
h = imrect;
position = round(wait(h));
close
xrange=position(1):position(1)+position(3);
yrange=position(2):position(2)+position(4);
reorientdataT1=reorientdataT1(yrange,xrange,:,:,:);

figure, imagesc(reorientdataT1(:,:,round(size(reorientdataT1,3)/2),1,1))


sl=':';
sl=5:17
dataT1=reorientdataT1(:,:,sl,:,:);


 imagesize=size(dataT1);
 imagesize=imagesize(1:end-1);
% new path
filepath_T1data=cell2mat(folderfiles{1}(1)); 
%% Fitting
       
       tic
       [fitparam]=T1FITTING_Pat_Aug2018(dataT1,TIs(:));
       toc
    
    
    
 



    

%% R^2 calculation
myfval=reshape(fitparam.fval,size(dataT1));

dimTI=find(numel(TIs)==size(dataT1));

ymean=1/numel(TIs)*sum(dataT1,dimTI);

repsize=ones(1,ndims(dataT1));
repsize(dimTI)=numel(TIs);

ymean_mat=repmat(ymean,repsize);

SStot=sum((dataT1-ymean_mat).^2,dimTI);

SSres=sum((myfval-dataT1).^2,dimTI);



R2=1-SSres./SStot;

R2adjusted=1-(numel(TIs)-1)/(numel(TIs)-3-1).*(1-R2);



R2adjusted_map=reshape(R2adjusted,imagesize(1),imagesize(2)*imagesize(3));

figure,imagesc(R2adjusted_map),caxis([0.5 1]),title('Adjusted R^2 map   (For Inversion Recovery fit)'), colorbar

%% mask out data that fit not a good fit.
myr2lim=0.8;
rsquarefilt=R2adjusted>myr2lim;
rsquarefilt=double(rsquarefilt);
rsquarefilt(rsquarefilt==0)=nan;


T1s=fitparam.RelaxTime.*rsquarefilt;
M0s=fitparam.M0.*rsquarefilt;
FlipAngles=fitparam.FlipAngle.*rsquarefilt;


%% reshape for plots

 
T1maps=reshape(T1s,imagesize(1),imagesize(2)*imagesize(3));
M0maps=reshape(M0s,imagesize(1),imagesize(2)*imagesize(3));
FAmaps=reshape(FlipAngles,imagesize(1),imagesize(2)*imagesize(3));
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
 
             
 fullsizeim=flip(rot90(double(data_cell{1}.img)),2); %open original image to substitute with t1maps values (xray effect)
  fullsizeim=fullsizeim*0.5; %fill the original image with t1maps values (xray effect)
  fullsizeim(yrange,xrange,sl)=tosave_T1maps;
  fullsizeim=rot90(flip(fullsizeim,2),-1); %Re-orient to initial orientation

 
  
 mynii.img=fullsizeim;
 
 
 %make sure to go one dimension less as we now have maps
 %niiheader.hdr.dime.dim(2:5) = [s1,s2,s3,s4];
 %niiheader.hdr.dime.pixdim(5) = 0;

 %make sure output is float
 %niiheader.hdr.dime.datatype=16;
 %niiheader.hdr.dime.bitpix=32;
 outfilename = [the_folder,slsh,the_filename,'___T1MAPS_test'];
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
T11=squeeze(T1maps(py,px,:,:));
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