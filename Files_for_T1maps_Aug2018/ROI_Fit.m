function [T1cal ,T1conf_low ,T1conf_high] = ROI_Fit(TIs, dataT1_ROI)
%ROI_Fit(TIS_REP,DATAT1_ROI_LINE)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : TIs
%      Y Output: dataT1_ROI
%  Output:
%      T1cal: T1 fit value
%      T1conf_low: lower bound T1 with 95% confidence
%       T1conf_high : Upper bound T1 with 95% confidence
%     
% 
TIs=TIs(:);
dataT1_ROI_line=reshape(dataT1_ROI,size(dataT1_ROI,1)*size(dataT1_ROI,2),size(dataT1_ROI,3));

normfact=nanmean(dataT1_ROI_line,2);
dataT1_ROI_line2=dataT1_ROI_line./(10*repmat(normfact,[1,size(dataT1_ROI_line,2)]));

np=size(dataT1_ROI_line2,1);

% TIs_rep=repmat(TIs,[size(dataT1_ROI_line,1),1]);
% 
% 
% %% Fit: 'untitled fit 1'.
% [xData, yData] = prepareCurveData( TIs_rep, dataT1_ROI_line2 );
% 
% % Set up fittype and options.
% ft = fittype( 'a*abs(1-(1+b)*exp(-x/c))', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% %opts.Robust = 'Bisquare';
% opts.Display = 'Off';
% opts.Lower = [0 0 200];
% opts.StartPoint = [0.5 1 0.8];
% opts.Upper = [1 1 5000];
% 
% % Fit model to data.
% [fitresultAV, gofAV] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% figure( 'Name', 'Averaged T1 fit' );
% h = plot( fitresultAV, xData, yData );
% legend( h, 'IR-TSE data vs TI ', 'T1 fit', 'Location', 'NorthEast' );
% % Label axes
% xlabel 'TI [ms]'
% ylabel 'Signal intensity [a.u.]'
% grid on


%% voxel by voxel

parfor ii=1:np
    
    
    
    ft = fittype( 'a*abs(1-(1+b)*exp(-x/c))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%opts.Robust = 'Bisquare';
opts.Display = 'Off';
opts.Lower = [0 0 200];
opts.StartPoint = [0.5 1 0.8];
opts.Upper = [1 1 5000];

    
    
    datay=dataT1_ROI_line2(ii,:);
    xData = TIs; 
    yData = datay(:);
    if ~any(yData)
        T1cal(ii)=nan;
        T1conf_high(ii)=nan;
        T1conf_low(ii)=nan;
    else
        [fitresult1, ~] = fit( xData, yData, ft, opts );
        coeff=coeffvalues(fitresult1);
        conf=confint(fitresult1);
        T1cal(ii)=coeff(3);
        T1conf_high(ii)=conf(2,3);
        T1conf_low(ii)=conf(1,3);
    end

  
  prog=(100*ii/np);
  if mod(prog,1)<0.01;
 disp([num2str(round(prog)),'%'])
  end
 
end
T1cal(T1cal<0)=nan;
T1cal(T1cal>4500)=nan;

T1conf_high(T1conf_high<100)=nan;
T1conf_high(T1conf_high>5000)=nan;

T1conf_low(T1conf_low<100)=nan;
T1conf_low(T1conf_low>5000)=nan;


T1cal=reshape(T1cal,[size(dataT1_ROI,1),size(dataT1_ROI,2)]);
 T1conf_high=reshape(T1conf_high,[size(dataT1_ROI,1),size(dataT1_ROI,2)]);
 T1conf_low=reshape(T1conf_low,[size(dataT1_ROI,1),size(dataT1_ROI,2)]);
 


