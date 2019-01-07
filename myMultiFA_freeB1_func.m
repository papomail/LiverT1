function [fitparam1, fitparam2, fitparam3, myfeval]= myMultiFA_freeB1_func(yData,xData,TR,options)
% fitparam1(=M0), fitparam2(=TR), myfeval
%   yData(=multiFA data in a voxel), xData(=FAs in radians), TR
% MaxIter=2*5400;
% 
% options.Display='notify';
% options.MaxFunEvals= '200*numberofvariables';
% options.MaxIter= MaxIter;
% options.TolFun= 1.000e-0;
% options.TolX= 1.000e-0;
% options.FunValCheck= 'off';
% options.OutputFcn=[];       
               

maxyy=max(yData);
 
init_val=[25*maxyy,1000,1];%
lowerbound=[5*maxyy,300,1];
upperbound=[500*maxyy,3000,1];


fhandle =@(x) norm( yData - (x(1).*sin(xData*x(3)).*(1-exp(-TR/x(2)))./(1-cos(xData*x(3)).*exp(-TR/x(2))))  );
%fhandle =@(x) norm( yData - x(1).*xData./(1+xData.^2.*x(2)/(2*TR))  ); %% Taylor expansion of the same function. However, it's not faster.

%[fitparam] = fminsearch(fhandle,init_val);
%[fitparam] = fminsearch(fhandle,init_val,options);

[fitparam] = fminsearchbnd(fhandle,init_val,lowerbound,upperbound,options);

fitparam1=fitparam(1);
fitparam2=fitparam(2);
fitparam3=fitparam(3);

myfeval=fitparam1.*sin(xData*fitparam3).*(1-exp(-TR/fitparam2))./(1-cos(xData*fitparam3).*exp(-TR/fitparam2));
end

