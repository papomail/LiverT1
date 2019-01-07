function [fitparam1, fitparam2, myfeval,fitparam3]= myhfunc(yData,xData,TR)
% fitparam1(=M0), fitparam2(=TR), myfeval
%   yData(=multiFA data in a voxel), xData(=FAs in radians), TR

options.Display='notify';
options.MaxFunEvals= '200*numberofvariables';
options.MaxIter= '200*numberofvariables';
options.TolFun= 1.000000000000000e-04;
options.TolX= 1.000000000000000e-04;
options.FunValCheck= 'off';
options.OutputFcn=[];       
               

maxyy=max(yData);
 
init_val=[100*maxyy,1000,-17000000];%
lowerbound=[maxyy,100,-5000*maxyy];
upperbound=[5000*maxyy,4000,5000*maxyy];


fhandle =@(x) norm( yData - (x(3)+x(1).*sin(xData).*(1-exp(-TR/x(2)))./(1-cos(xData).*exp(-TR/x(2))))  );
%fhandle =@(x) norm( yData - x(1).*xData./(1+xData.^2.*x(2)/(2*TR))  ); %% Taylor expansion of the same function. However, it's not faster.

%[fitparam] = fminsearch(fhandle,init_val);
%[fitparam] = fminsearch(fhandle,init_val,options);

[fitparam] = fminsearchbnd(fhandle,init_val,lowerbound,upperbound,options);

fitparam1=fitparam(1);
fitparam2=fitparam(2);
fitparam3=fitparam(3);
myfeval=fitparam3+fitparam1.*sin(xData).*(1-exp(-TR/fitparam2))./(1-cos(xData).*exp(-TR/fitparam2));
end

