function [ R2adjusted ] = R2calc( dataT1,Xvals,myfval )

imagesize=size(dataT1);
 imagesize=imagesize(1:end-1);
% Rsquare and AdjustedRsquare calculation
%

numX=numel(Xvals);

%myfval=reshape(fitparam.fval,size(dataT1));

dimTI=find(numX==size(dataT1));

ymean=1/numX*sum(dataT1,dimTI);

repsize=ones(1,ndims(dataT1));
repsize(dimTI)=numX;

ymean_mat=repmat(ymean,repsize);

SStot=sum((dataT1-ymean_mat).^2,dimTI);

SSres=sum((myfval-dataT1).^2,dimTI);



R2=1-SSres./SStot;

R2adjusted=1-(numX-1)/(numX-3-1).*(1-R2);



R2adjusted_map=reshape(R2adjusted,imagesize(1),imagesize(2)*imagesize(3));

if Xvals(end)<1000
    titletext='Adjusted R^2 map   (For multi FA fit)';
elseif Xvals(end)>1000
    titletext='Adjusted R^2 map   (For Inversion Recovery fit)';
    
end
figure,imagesc(R2adjusted_map),caxis([0.5 1]),title(titletext), colorbar

end

