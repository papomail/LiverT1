%dirnames{1}=the_folder;

%% Do the fitting


%% fitting Function

function [output]=T1FITTING_MULTIFA_Pat_Aug2018(data,xvals)
%function [fitparam,ci]=fit_exp_decay(data,xvals)
xdim=find(numel(xvals)==size(data));
   
if xdim==3 && ndims(data)==4
            data=permute(data,[1 2 4 3]);
            ydata=reshape(data,[size(data,1)*size(data,2)*size(data,3),size(data,4)]);

elseif xdim==4 && ndims(data)==5
            data=permute(data,[1 2 3 5 4]);
            ydata=reshape(data,[size(data,1)*size(data,2)*size(data,3)*size(data,4),size(data,5)]);

elseif xdim==5 && ndims(data)==5            
            ydata=reshape(data,[size(data,1)*size(data,2)*size(data,3)*size(data,4),size(data,5)]);
                        
elseif xdim==4 && ndims(data)==4           
            ydata=reshape(data,[size(data,1)*size(data,2)*size(data,3),size(data,4)]);
else 
    CHECK_YOUR_DIMENSIONS
end
Tmean=mean(xvals);
tmin=min(xvals);
options=optimset('Display','off');

M0=zeros(size(ydata,1),1);
aa=M0;

RelaxTime=zeros(size(ydata,1),1);
        parfor ii=1:size(ydata,1)
		z = squeeze(ydata(ii,:));
		z=z(:);
        %init_val=[2*max(z),Tmean,pi];% [x(1),x(2),x(3)
		%fhandle = @(x) norm(z - abs( x(1)* (1+ (cos(x(3))-1)*exp(-xvals./x(2)) )) );
        %fhandle = @(x) norm(z - abs(x(1)*(1    -x(3)*exp(-xvals./x(2)))));
       
        %init_val=[0,Tmean,-max(z),1.2*max(z)];% [x(1),x(2),x(3),x(4)
        %fhandle = @(x) norm(z -   ( x(1)+abs(x(3)*exp(-xvals./x(2))+x(4)*(1-exp(-xvals./x(2)))) )   );

        
        init_val=[1.3*max(z),1,Tmean];% [x(1),x(2),x(3),x(4)
       
        fhandle = @(x) norm(z -   ( x(1)*abs(1-(1+x(2))*exp(-xvals/x(3))) )   );
      
		[fitparam] = fminsearch(fhandle,init_val,options);
        
        
        
       
%         
%         aa(ii) = fitparam(1);
		RelaxTime(ii) = fitparam(3);
        FlipAngle(ii) = fitparam(2);
         M0(ii)=fitparam(1);
         
         myfval(ii,:)= fitparam(1)*abs(1-(1+fitparam(2))*exp(-xvals/fitparam(3)));
        end
        
        if xdim==3
              M0=reshape(M0,[size(data,1),size(data,2),size(data,3)]);
              RelaxTime=reshape(RelaxTime,[size(data,1),size(data,2),size(data,3)]);
              FlipAngle=reshape(FlipAngle,[size(data,1),size(data,2),size(data,3)]);
              myfval=reshape(myfval,[size(data,1),size(data,2),size(data,3)]);

        elseif xdim==4 && ndims(data)==5
              M0=reshape(M0,[size(data,1),size(data,2),size(data,3),size(data,4)]);
              RelaxTime=reshape(RelaxTime,[size(data,1),size(data,2),size(data,3),size(data,4)]);
              FlipAngle=reshape(FlipAngle,[size(data,1),size(data,2),size(data,3),size(data,4)]);
              myfval=reshape(myfval,[size(data,1),size(data,2),size(data,3),size(data,4),size(data,5)]);

        elseif  xdim==4 && ndims(data)==4  
                 %       aa=reshape(aa,[size(data,1),size(data,2),size(data,3)]);
            myfval=reshape(myfval,[size(data,1),size(data,2),size(data,3),size(data,4)]);
            M0=reshape(M0,[size(data,1),size(data,2),size(data,3)]);
              RelaxTime=reshape(RelaxTime,[size(data,1),size(data,2),size(data,3)]);
              FlipAngle=reshape(FlipAngle,[size(data,1),size(data,2),size(data,3)]);
              
       
        
        else xdim==5 && ndims(data)==5            
             % aa=reshape(aa,[size(data,1),size(data,2),size(data,3),size(data,4)]);
            myfval=reshape(myfval,[size(data,1),size(data,2),size(data,3),size(data,4),size(data,5)]);
            M0=reshape(M0,[size(data,1),size(data,2),size(data,3),size(data,4)]);
              RelaxTime=reshape(RelaxTime,[size(data,1),size(data,2),size(data,3),size(data,4)]);
              FlipAngle=reshape(FlipAngle,[size(data,1),size(data,2),size(data,3),size(data,4)]);         
        
        end
        %output.aa=aa;
        output.M0=M0;
        output.RelaxTime=RelaxTime;
        output.FlipAngle=FlipAngle;
        
       
        
        output.fval=myfval;
%         
        
        
        
        
        
        
        
        
end
% 

