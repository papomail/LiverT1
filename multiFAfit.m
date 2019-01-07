function [output] = multiFAfit(dataT1,FAs,TR,varargin)
dataT1=double(dataT1);
FAs=double(FAs);


%clear T1cal T1conf_high T1conf_low
xData=pi/180*FAs(:);
%coder.extrinsic('myhfunc')

mysize=size(dataT1);
linesize=[prod(mysize(1:end-1)) mysize(end)];

dataT1line=reshape(dataT1,linesize);


if nargin<4
    warning('B1 MAP not found. Will continue without it but the T1 maps WONT be corrected for B1 inhomogeneities')
    dataB1line=100*ones(linesize(1),1);

else
   dataB1=varargin{1};
   dataB1line=reshape(dataB1,linesize(1),1);
   
end
 







RelaxTime=zeros(linesize(1:end-1),1); 
M0=zeros(linesize(1:end-1),1);
%Aconstant=M0;
nv=linesize(1);
myfval=zeros(linesize);
FA_corrected=myfval;


MaxIter=5000;

options.Display='notify';
options.MaxFunEvals= '200*numberofvariables';
    options.MaxFunEvals= 5000;

options.MaxIter= MaxIter;
options.TolFun= 1.000e-0;
options.TolX= 1.000e-0;
options.FunValCheck= 'off';
options.OutputFcn=[];   

parfor ii=1:nv
   
 b1=dataB1line(ii);       
  xData_cor=xData*b1/100;

yData=dataT1line(ii,:);
yData=yData(:);


%fitparam=0;

%[fitparam1, fitparam2, feval] = myMultiFAfunc(yData,xData_cor,TR,options);
%[fitparam1, fitparam2, feval,fitparam3] = myMultiFAfunc_withConstant(yData,xData_cor,TR,options);
[fitparam1, fitparam2, fitparam3, feval] = myMultiFA_freeB1_func(yData,xData_cor,TR,options);



%  init_val=[100*maxyy,1000];% [x(1),x(2),x(3),x(4)
%     
%      
% 
% fhandle =@(x) norm( yData - x(1).*sin(xData).*(1-exp(-4/x(2)))./(1-cos(xData).*exp(-4/x(2)))  );
% [fitparam] = fminsearch(fhandle,init_val,options);
%          %ft = fittype( 'a*x*(1+4/b/2)/(1+4/b/2+x^2/2/(4/b))', 'independent', 'x', 'dependent', 'y' );

        M0(ii)=fitparam1;
        RelaxTime(ii) = fitparam2;
        myfval(ii,:) = feval;
        FA_corrected(ii,:) = xData_cor*180/pi;

     %   Aconstant(ii)=fitparam3;
     B1map(ii)=fitparam3;

end
   
    
RelaxTime=reshape(RelaxTime,mysize(1:end-1));
M0=reshape(M0,mysize(1:end-1));
myfval=reshape(myfval,mysize);
FA_corrected=reshape(FA_corrected,mysize);
%Aconstant=reshape(Aconstant,mysize(1:end-1));
B1map=reshape(B1map,mysize(1:end-1));

 output.M0=M0;
 output.RelaxTime=RelaxTime;
 output.fval=myfval;
 output.FA_b1corrected=FA_corrected;
% output.Aconstant=Aconstant;
output.B1map=B1map;
end
