
function PlotT1Fits(T1maps,FAmaps,M0maps, use_b1_in_fit)
global T1maps FAmaps M0maps use_b1_in_fit
mouseMove
    function mouseMove
        %  global  The_CEST_data
        global T1maps FAmaps M0maps use_b1_in_fit
        
        Ysize=size(M0maps,1);
        Xsize=size(M0maps,2);
        
        
        figure(1)
        C = get (gca, 'CurrentPoint');
        %title(gca, ['(X,Y,Z) = (', num2str(round(C(1,1))), ', ',num2str(round(C(1,2))),', ',sprintf('%.0f%',slicez),')']);
        
        
        
        
        xx=round(C(1,1));
        yy=round(C(1,2));
        
        if xx>0 && xx<Xsize+1 && yy>0 && yy<Ysize+1
            
            %Plot selected voxels  (only works for IR for the moment)
            px=round(xx);
            py=round(yy);
            
            %aa1=squeeze(aamaps(py,px,:,:));
            M01=squeeze(M0maps(py,px,:,:));
            if exist('LookLocker','var') && strcmp(LookLocker,'yes')
                T11=squeeze(T1LLfit(py,px,:,:));
            else
                T11=squeeze(T1maps(py,px,:,:));
            end
            
            reshaped_dataT1=reshape(dataT1,size(dataT1,1),size(dataT1,2)*size(dataT1,3),size(dataT1,4),size(dataT1,5));
            dataT11=squeeze(reshaped_dataT1(py,px,:,:));
            
            
            switch fit_type
                case 'inversion_recovery'
                    
                    figure,
                    plot(TIs,squeeze(dataT11(:,1)),'o')
                    hold on
                    %plot(myxx,aa1(1)+abs(M01(1)*exp(-myxx/T11(1))+Minf1(1)*(1-exp(-myxx/T11(1) ))))
                    
                    if use_b1_in_fit
                        B1map=reshape(dataB1,size(dataB1,1),size(dataB1,2)*size(dataB1,3));
                        Minf1=-cos(pi*B1map(py,px,:,:)/100);
                    else
                        Minf1=squeeze(FAmaps(py,px,:,:));
                    end
                    
                    myxx=logspace(0,3.8,300);
                    plot(myxx,M01*abs(1-(1+Minf1)*exp(-myxx/T11)))
                    xlabel('Inversion delay (ms)')
                    ylabel('Signal Intensity')
                case 'multi_fa'
                    
                    reshaped_FAcor=reshape(fitparam.FA_b1corrected,size(dataT1,1),size(dataT1,2)*size(dataT1,3),size(dataT1,4),size(dataT1,5));
                    FAcor1=squeeze(reshaped_FAcor(py,px,:,:));
                    
                    myxx=linspace(0,20,300);
                    figure,
                    plot(FAcor1,squeeze(dataT11(:,1)),'o')
                    hold on
                    
                    if exist('Aconstantmap','var')
                        Aconstantmap1=squeeze(Aconstantmap(py,px,:,:));
                    else
                        Aconstantmap1=0;
                    end
                    plot(myxx,Aconstantmap1+M01*sin(myxx*pi/180).*(1-exp(-TR/T11))./(1-cos(myxx*pi/180).*exp(-TR/T11)))
                    xlabel('Efective FlipAngle [degrees]')
                    ylabel('Signal Intensity')
                    
                    
            end
            
        end
        
        
    end
end



