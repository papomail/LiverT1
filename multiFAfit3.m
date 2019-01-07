function [RelaxTime, M0] = multiFAfit3(dataT1,FAs,TR)
dataT1=double(dataT1);
FAs=double(FAs);


%clear T1cal T1conf_high T1conf_low
xData=pi/180*FAs(:);
coder.extrinsic('myhfunc')



mysize=size(dataT1);
linesize=[prod(mysize(1:end-1)) mysize(end)];

dataT1line=reshape(dataT1,linesize);

RelaxTime=zeros(linesize(1:end-1),1); 
M0=zeros(linesize(1:end-1),1);
nv=linesize(1);

parfor ii=1:nv
   
        
        

yData=dataT1line(ii,:);
yData=yData(:);


%fitparam=0;

[fitparam1, fitparam2] = myhfunc(yData,xData,TR);
%  init_val=[100*maxyy,1000];% [x(1),x(2),x(3),x(4)
%     
%      
% 
% fhandle =@(x) norm( yData - x(1).*sin(xData).*(1-exp(-4/x(2)))./(1-cos(xData).*exp(-4/x(2)))  );
% [fitparam] = fminsearch(fhandle,init_val,options);
%          %ft = fittype( 'a*x*(1+4/b/2)/(1+4/b/2+x^2/2/(4/b))', 'independent', 'x', 'dependent', 'y' );

        RelaxTime(ii) = fitparam2;
         M0(ii)=fitparam1;



end
   
    
RelaxTime=reshape(RelaxTime,mysize(1:end-1));
M0=reshape(M0,mysize(1:end-1));

end
