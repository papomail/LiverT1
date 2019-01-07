function [RelaxTime, M0] = multiFAfit2(dataT1,FAs,TR)
dataT1=double(dataT1);
FAs=double(FAs);

RelaxTime=zeros([size(dataT1,1),size(dataT1,2)]); %%% TEMP!!
M0=zeros([size(dataT1,1),size(dataT1,2)]);
%clear T1cal T1conf_high T1conf_low
xData=pi/180*FAs(:);
coder.extrinsic('myhfunc')



mysize=size(dataT1)
prod(mysize(1:end-1))

for ii=1:size(dataT1,1)
    for jj=1:size(dataT1,2)
        
        

y=dataT1(ii,jj,2,1,:);
yData=y(:);


%fitparam=0;

[fitparam1, fitparam2] = myhfunc(yData,xData,TR);
%  init_val=[100*maxyy,1000];% [x(1),x(2),x(3),x(4)
%     
%      
% 
% fhandle =@(x) norm( yData - x(1).*sin(xData).*(1-exp(-4/x(2)))./(1-cos(xData).*exp(-4/x(2)))  );
% [fitparam] = fminsearch(fhandle,init_val,options);
%          %ft = fittype( 'a*x*(1+4/b/2)/(1+4/b/2+x^2/2/(4/b))', 'independent', 'x', 'dependent', 'y' );

        RelaxTime(ii,jj) = fitparam2;
         M0(ii,jj)=fitparam1;



    end
   
end
