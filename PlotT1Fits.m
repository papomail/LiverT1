
function PlotT1Fits(T1maps,FAmaps,M0maps,TIs,dataT1,dataB1, use_b1_in_fit,fit_type,R2adjusted,fitparam,TR)
fig1=figure;
%subplot(2,1,1)
h=subaxis(2,2,1,1, 'm', 0.07);
imagesc(T1maps),title('T1 map (ms). Pan through to see the fits'),
colorbar
caxis([0 4000])
colormap(h,jet)
set(gca,'XTick','')
set(gca,'YTick','')
set(gca,'FontSize',12);


med_M0=nanmedian(M0maps(:));
Ysize=size(M0maps,1);
Xsize=size(M0maps,2);

r2amap=reshape(R2adjusted,[size(R2adjusted,1),size(R2adjusted,2)*size(R2adjusted,3)]);
h2=subaxis(2,2,1,2, 'margin', 0.07);
imagesc(r2amap),caxis([0.6 1])
colormap(h2,parula)
colorbar
xlabel('Adjusted R^2')
set(gca,'XTick','')
set(gca,'YTick','')
set(gca,'FontSize',12);


if isfield(fitparam,'FA_b1corrected')
reshaped_FAcor=reshape(fitparam.FA_b1corrected,size(dataT1,1),size(dataT1,2)*size(dataT1,3),size(dataT1,4),size(dataT1,5));
end

B1map=reshape(dataB1,size(dataB1,1),size(dataB1,2)*size(dataB1,3),size(dataB1,4),size(dataB1,5));

myxx=linspace(0,40,100);


mouseMove;
set (fig1, 'WindowButtonMotionFcn', @mouseMove,'pointer','crosshair')




    function mouseMove(~,~)
        
        C = get (h, 'CurrentPoint');
        C2 = get (h2, 'CurrentPoint');
        
        %title(gca, ['(X,Y,Z) = (', num2str(round(C(1,1))), ', ',num2str(round(C(1,2))),', ',sprintf('%.0f%',slicez),')']);
        
        xx=round(C(1,1));
        yy=round(C(1,2));
        xx2=round(C2(1,1));
        yy2=round(C2(1,2));
        if (xx>0 && xx<Xsize+1 && yy>0 && yy<Ysize+1)
            
            px=round(xx);
            py=round(yy);
            do_theplots
        elseif (xx2>0 && xx2<Xsize+1 && yy2>0 && yy2<Ysize+1)
            px=round(xx2);
            py=round(yy2);
            do_theplots
        end
        
        function do_theplots(~ ,~)
            %aa1=squeeze(aamaps(py,px,:,:));
            M01=squeeze(M0maps(py,px,:,:));
            if exist('LookLocker','var') && strcmp(LookLocker,'yes')
                T11=squeeze(T1LLfit(py,px,:,:));
            else
                T11=squeeze(T1maps(py,px,:,:));
            end
            
            reshaped_dataT1=reshape(dataT1,size(dataT1,1),size(dataT1,2)*size(dataT1,3),size(dataT1,4),size(dataT1,5));
            dataT11=squeeze(reshaped_dataT1(py,px,:,:));
            
            subaxis(2,2,2,1,1.1,1.95,'pl', 0.05);
            switch fit_type
                case 'inversion_recovery'
                    yydata=squeeze(dataT11(:,1));
                    
                    %subplot(2,1,2)
                    %subaxis(2, 1, 2, 'margin', 0.08);
                    
                    
                    %plot(myxx,aa1(1)+abs(M01(1)*exp(-myxx/T11(1))+Minf1(1)*(1-exp(-myxx/T11(1) ))))
                    
                    if use_b1_in_fit
                        Minf1=-cos(pi*B1map(py,px,:,:)/100); 

                    else
                        Minf1=squeeze(FAmaps(py,px,:,:));
                    end
                    
                    myxx=logspace(0,3.8,300);
                    plot(myxx,M01*abs(1-(1+Minf1)*exp(-myxx/T11)),'LineWidth',1.5)
                    hold on
                    plot(TIs,yydata,'o','LineWidth',2)
                    xlim([0 3000]);
                    %ylim([0 1.5*med_M0]);
                    
                    
                    %title('T1 map [ms]. Pan through to see the fits')
                    xlabel('Inversion delay (ms)')
                    ylabel('Signal Intensity')
                    %legend('Data','Fit')
                    mytext=['T1=',mat2str(round(T11)),' ms'];
                    text(0.5,0.8,mytext,'units','normalized','FontSize',14)
                    set(gca,'FontSize',12);
                    
                    hold off;
                    
                case 'multi_fa'
                    
                    if exist('Aconstantmap','var')
                        Aconstantmap1=squeeze(Aconstantmap(py,px,:,:));
                    else
                        Aconstantmap1=0;
                    end
                    FAcor1=squeeze(reshaped_FAcor(py,px,:,:));
                    
                    cla
                    
                    fplot(@(x) Aconstantmap1+M01*sin(x*pi/180).*(1-exp(-TR/T11))./(1-cos(x*pi/180).*exp(-TR/T11)), [0 20])
                    hold on
                    plot(FAcor1,squeeze(dataT11(:,1)),'o')
                    %plot(myxx,Aconstantmap1+M01*sin(myxx*pi/180).*(1-exp(-TR/T11))./(1-cos(myxx*pi/180).*exp(-TR/T11)))
                    xlabel('Efective FlipAngle [degrees]')
                    ylabel('Signal Intensity')
                    mytext=['T1=',mat2str(round(T11)),' ms'];
                    text(0.5,0.8,mytext,'units','normalized','FontSize',14)
                    set(gca,'FontSize',12);
                   
            end
        end
    end




end



