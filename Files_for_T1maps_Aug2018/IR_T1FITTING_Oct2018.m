%dirnames{1}=the_folder;


%% fitting Function

function [output]=IR_T1FITTING_Oct2018(dataT1,xvals,varargin)
%function [fitparam,ci]=fit_exp_decay(data,xvals)
xdim=find(numel(xvals)==size(dataT1));





if nargin<3
    warning('B1 MAP not found. Will continue without it but the T1 maps WONT be corrected for B1 inhomogeneities')
    %dataB1line=ones(linesize(1),1);
    
else
    
mysize=size(dataT1);
linesize=[prod(mysize(1:end-1)) mysize(end)];

    dataB1=varargin{1};
    dataB1line=reshape(dataB1,linesize(1),1);
    
end






if xdim==3 && ndims(dataT1)==4
    dataT1=permute(dataT1,[1 2 4 3]);
    ydata=reshape(dataT1,[size(dataT1,1)*size(dataT1,2)*size(dataT1,3),size(dataT1,4)]);
    
elseif xdim==4 && ndims(dataT1)==5
    dataT1=permute(dataT1,[1 2 3 5 4]);
    ydata=reshape(dataT1,[size(dataT1,1)*size(dataT1,2)*size(dataT1,3)*size(dataT1,4),size(dataT1,5)]);
    
elseif xdim==5 && ndims(dataT1)==5
    ydata=reshape(dataT1,[size(dataT1,1)*size(dataT1,2)*size(dataT1,3)*size(dataT1,4),size(dataT1,5)]);
    
elseif xdim==4 && ndims(dataT1)==4
    ydata=reshape(dataT1,[size(dataT1,1)*size(dataT1,2)*size(dataT1,3),size(dataT1,4)]);
else
    CHECK_YOUR_DIMENSIONS
end
%Tmean=mean(xvals);
%tmin=min(xvals);
options=optimset('Display','off');

M0=zeros(size(ydata,1),1);
%aa=M0;

RelaxTime=zeros(size(ydata,1),1);


ydata=double(ydata);

if exist('dataB1line','var') %% FIT with B1 
    
    
MaxIter=5000; %Only two parameters to fit (M0, T1)
options.Display='notify';
options.MaxFunEvals= '200*numberofvariables';
options.MaxIter= MaxIter;
options.TolFun= 1.00e-0;
options.TolX= 1.00e-0;
options.FunValCheck= 'off';
options.OutputFcn=[];       
    
    parfor ii=1:size(ydata,1)
        z = squeeze(ydata(ii,:));
        z=z(:);
        b1=pi*dataB1line(ii)/100;
        maxyy=max(z);
        %init_val=[2*max(z),Tmean,pi];% [x(1),x(2),x(3)
        fhandle = @(x) norm(z - abs( x(1)* (1+ (cos(b1)-1)*exp(-xvals./x(2)) )) );
        %fhandle = @(x) norm(z - abs(x(1)*(1    -x(3)*exp(-xvals./x(2)))));
        
        %init_val=[0,Tmean,-max(z),1.2*max(z)];% [x(1),x(2),x(3),x(4)
        %fhandle = @(x) norm(z -   ( x(1)+abs(x(3)*exp(-xvals./x(2))+x(4)*(1-exp(-xvals./x(2)))) )   );
        [~,per]=sort(z);
        nullti=nanmean(xvals(per(1:3)));
        initT1=nullti/log(2);
        init_val=[1.3*maxyy,initT1];% [x(1),x(2),x(3),x(4)
        lowerbound=[maxyy,100];
        upperbound=[50*maxyy,4000];
        %fhandle = @(x) norm(z -   ( x(1)*abs(1-(1+x(2))*exp(-xvals/x(3))) )   );
        %[fitparam] = fminsearch(fhandle,init_val,options);
        [fitparam] = fminsearchbnd(fhandle,init_val,lowerbound,upperbound,options);
        
        %         aa(ii) = fitparam(1);
        RelaxTime(ii) = fitparam(2);
        M0(ii)=fitparam(1);
        
        myfval(ii,:)= fitparam(1)*abs(1+ (cos(b1)-1)*exp(-xvals/fitparam(2)));
    end
    
    
    
    if xdim==3
        M0=reshape(M0,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        RelaxTime=reshape(RelaxTime,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        myfval=reshape(myfval,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        
    elseif xdim==4 && ndims(dataT1)==5
        M0=reshape(M0,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        RelaxTime=reshape(RelaxTime,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        myfval=reshape(myfval,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4),size(dataT1,5)]);
        
    elseif  xdim==4 && ndims(dataT1)==4
        %       aa=reshape(aa,[size(data,1),size(data,2),size(data,3)]);
        myfval=reshape(myfval,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        M0=reshape(M0,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        RelaxTime=reshape(RelaxTime,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        
        
        
    else xdim==5 && ndims(dataT1)==5
        % aa=reshape(aa,[size(data,1),size(data,2),size(data,3),size(data,4)]);
        myfval=reshape(myfval,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4),size(dataT1,5)]);
        M0=reshape(M0,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        RelaxTime=reshape(RelaxTime,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        
    end
    %output.aa=aa;
    output.M0=M0;
    output.RelaxTime=RelaxTime;
    output.fval=myfval;
    
    
    
    
else %% if no b1 is used in the fit
    
    
MaxIter=8000; %3 parameters to fit (M0,FA, T1)
options.Display='notify';
options.MaxFunEvals= '200*numberofvariables';
options.MaxIter= MaxIter;
options.TolFun= 1.00e-0;
options.TolX= 1.000e-0;
options.FunValCheck= 'off';
options.OutputFcn=[];       
    parfor ii=1:size(ydata,1)
        z = squeeze(ydata(ii,:));
        z=z(:);
        maxyy=max(z);
        %init_val=[2*max(z),Tmean,pi];% [x(1),x(2),x(3)
        %fhandle = @(x) norm(z - abs( x(1)* (1+ (cos(x(3))-1)*exp(-xvals./x(2)) )) );
        %fhandle = @(x) norm(z - abs(x(1)*(1    -x(3)*exp(-xvals./x(2)))));
        
        %init_val=[0,Tmean,-max(z),1.2*max(z)];% [x(1),x(2),x(3),x(4)
        %fhandle = @(x) norm(z -   ( x(1)+abs(x(3)*exp(-xvals./x(2))+x(4)*(1-exp(-xvals./x(2)))) )   );
        [~,per]=sort(z);
        nullti=nanmean(xvals(per(1:3)));
        initT1=nullti/log(2);
        init_val=[1.3*maxyy,1,initT1];% [x(1),x(2),x(3),x(4)

        lowerbound=[maxyy,0,100];
        upperbound=[50*maxyy,1,4000];
        
        fhandle = @(x) norm(z -   ( x(1)*abs(1-(1+x(2))*exp(-xvals/x(3))) )   );
        [fitparam] = fminsearchbnd(fhandle,init_val,lowerbound,upperbound,options);
        
        
        %         aa(ii) = fitparam(1);
        RelaxTime(ii) = fitparam(3);
        FlipAngle(ii) = fitparam(2);
        M0(ii)=fitparam(1);
        
        myfval(ii,:)= fitparam(1)*abs(1-(1+fitparam(2))*exp(-xvals/fitparam(3)));
    end
    
    
    
    if xdim==3
        M0=reshape(M0,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        RelaxTime=reshape(RelaxTime,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        FlipAngle=reshape(FlipAngle,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        myfval=reshape(myfval,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        
    elseif xdim==4 && ndims(dataT1)==5
        M0=reshape(M0,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        RelaxTime=reshape(RelaxTime,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        FlipAngle=reshape(FlipAngle,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        myfval=reshape(myfval,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4),size(dataT1,5)]);
        
    elseif  xdim==4 && ndims(dataT1)==4
        %       aa=reshape(aa,[size(data,1),size(data,2),size(data,3)]);
        myfval=reshape(myfval,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        M0=reshape(M0,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        RelaxTime=reshape(RelaxTime,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        FlipAngle=reshape(FlipAngle,[size(dataT1,1),size(dataT1,2),size(dataT1,3)]);
        
        
        
    else xdim==5 && ndims(dataT1)==5
        % aa=reshape(aa,[size(data,1),size(data,2),size(data,3),size(data,4)]);
        myfval=reshape(myfval,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4),size(dataT1,5)]);
        M0=reshape(M0,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        RelaxTime=reshape(RelaxTime,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        FlipAngle=reshape(FlipAngle,[size(dataT1,1),size(dataT1,2),size(dataT1,3),size(dataT1,4)]);
        
    end
    %output.aa=aa;
    output.M0=M0;
    output.RelaxTime=RelaxTime;
    output.FlipAngle=FlipAngle;
    output.fval=myfval;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
end
%

