fclose all
clear all


%% Options
SHOW_IMAGES=0


%% Set original path
[the_folder] = pathsetting('Desktop/REVITA2');
%%  Study Folder Selection
message = 'Select the study folder. (Or individual cases) ';
folders=uipickfiles('FilterSpec',the_folder,'Prompt',message);

filter='mFAfit_T1MAPS_withB1.mat'
files = FilterFiles(folders,filter);
ncases=numel(files);

%% go through the cases
R2mean=zeros(1,ncases);
R2std=R2mean;
T1mean=R2mean;
T1std=R2mean;

for ii=1:ncases
    load(char(files(ii)));
    ss=size(dataT1);
    sm=ss(1:3);
    middlez=round(sm(3)/2);
    
     position=[round(0.2*ss(1))   round(0.35*ss(2))   round(ss(1)/4)   round(ss(2)/4)];
    xrange=position(1):position(1)+position(3);
    yrange=position(2):position(2)+position(4);    
    zrange=middlez-3:middlez+3;
    
    
    
    %% Liver parenchima --> take voxels within  { mean(T1) +- std(T1) }  (assumes the mayority of the ROI is parenchima)
    T1=reshape(T1maps,sm);
    sT1=(T1(yrange,xrange,zrange));
    t1mean=nanmean(sT1(:));
    t1std=nanstd(sT1(:));
    sT1g=sT1;
    sT1g(sT1g<(t1mean-t1std))=nan; 
    sT1g(sT1g>(t1mean+t1std))=nan;

    
        %% Take best fitted voxels  --> take voxels with Adjusted R squared  { mean(R2)-std(R2) and above }  

    sR2=R2adjusted(yrange,xrange,zrange);
    r2mean=nanmean(sR2(:));
    r2std=nanstd(sR2(:));
    sR2g=sR2;
    sR2g(sR2<(r2mean-r2std))=nan;
       R2mask= ~isnan(sR2g);
    

    %% Multiply R2good with T1good
    t1good=R2mask.*sT1g;
    
    
    %% maths
   T1mean(ii)=nanmean(t1good(:));
   T1std(ii)=nanstd(t1good(:));
   
   R2mean(ii)=nanmean(sR2g(:));
   R2std(ii)=nanstd(sR2g(:));
   
   
   
    %% show images
    if SHOW_IMAGES
    ca=[400 3000];
    figure,imagesc(T1(:,:,round(sm(3)/2))),colormap(jet),caxis(ca);title('T1')
     %figure,imagesc(sT1(:,:,4)),colormap(jet),caxis(ca);title('T1 in ROI')
     %figure,imagesc(sT1g(:,:,4)),colormap(jet),caxis(ca);title('Parenchyma T1')
     %figure,imagesc(sR2(:,:,4)),colormap(jet),caxis([.5 1]);title('aR^2')
     % figure,imagesc(sR2g(:,:,4)),colormap(jet),caxis([.5 1]);title('well fitted aR^2')

      figure,imagesc(t1good(:,:,4)),colormap(jet),caxis([.5 1]);caxis(ca);title('T1 parenchyma good  [ms]'),colorbar
    end
    
   
end



%% Plots
x = 1:ncases;
figure, 
bar(x,T1mean,'k')                

hold on

er = errorbar(x,T1mean,T1std*0,T1std/2);    
er.Color = 'k';                            
er.LineStyle = 'none';  
er.LineWidth = 2; 
title('Mean T1')
 


figure
bar(x,R2mean,'r')      
hold on
er = errorbar(x,R2mean,R2std*0,R2std/2);    
er.Color = 'r';                            
er.LineStyle = 'none';  
er.LineWidth = 2; 
 title('Mean aR^2')         
