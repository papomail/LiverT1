function showFatWater(II,varargin)
%% Show Fat and Water images. Usage: showFatWater(II,varargin), where II is the 4D(3D_fat,3D_water) data and varargin defines the zslice to show. If varargin is empty, then the midle zslize is
% shown
sl=round(size(II,3)/2);
FAscan=2;


showALLFAs=1
if showALLFAs
    FAscan=1:size(II,5);
else
    FAscan=1;
end

for ii=1:numel(FAscan)
    
    if ~isempty(varargin)
        sl=varargin{1};
    end
    %
    
    
    if size(II,4)==4
        
        figure,
        subplot(2,2,1),imagesc(II(:,:,sl,1,ii)),
        title('Fat only')
        set(gca,'XTick','')
        set(gca,'YTick','')
        subplot(2,2,2),imagesc(II(:,:,sl,2,ii)),
        title('Water only')
        set(gca,'XTick','')
        set(gca,'YTick','')
        subplot(2,2,3),imagesc(II(:,:,sl,3,ii)),
        title('Out of phase')
        set(gca,'XTick','')
        set(gca,'YTick','')
        subplot(2,2,4),imagesc(II(:,:,sl,4,ii)),
        title('In phase')
        set(gca,'XTick','')
        set(gca,'YTick','')
        colormap(gray)
        
    else
        
        % II=II(:,:,:,[1,4],:);  %% take FatOnly and WaterOnly images and get rid of the rest.
        
        figure,
        subplot(1,2,1),imagesc(II(:,:,sl,1,ii)),
        title('Water only')
        set(gca,'XTick','')
        set(gca,'YTick','')
        
        subplot(1,2,2),imagesc(II(:,:,sl,2,ii)),
        title('Fat only')
        set(gca,'XTick','')
        set(gca,'YTick','')
        colormap(gray)
    end
end

end
