%Renders an overlay of a fit parameter (ccv) over the photons array
%(photons). A handle to the graphics object is returned.

%The photons array will be scaled according to photonsMax
    %photonsMax = -1 -> scale LUTs to the maximum for the photon array
    %photonsMax = any positive integer -> scale photons to this as maxium
    %for consistent display of images in a set.
    
%The dimensions of photons and ccv must match.

%The color coded value (lifetime parameter) will be rendered in the ccvColor colormap across
%ccvRange, which should be specified as [limLow limHigh]
    
%The figure will display the title figTitle so that it is clear which image
%it corresponds to (especially important when it is called from
%renderImages or within the pxwiseAnalysis).

function [figHandle,f2,fImgOnly] = renderOverlay(photons,ccv,photonsMax,ccvRange,ccvColor,figTitle)

bright = 1; %if <1, it enhances photon brightness a bit to make the image clearer
cmap = colormap(ccvColor);
cmapFlip = flipud(cmap);

%size of input images
pSize = size(photons);
ccvSize = size(ccv);
if(size(ccv,3) ~= 1)
    error('Attempting to render overlay of a multi-D variable. Did you mean tm?');
end
if(~isequal(pSize,ccvSize))
    error('Attempt to render overlays of unequal size');
end

%lifetime - grayscale overlay
figHandle = figure('Name',figTitle);
axList = gca;
set(axList(1),'Visible','off');
hold on;
axList(2) = axes;
set(axList(2),'Visible','off');
hold on;
title(figTitle,'Visible','on','Interpreter','none');

%create a transparency map
if photonsMax == -1
    maxVal = max(photons,[],'all');
elseif photonsMax > 1
    maxVal = photonsMax;
else
    error('Parameter photonsMax should be either a positive integer or -1 (to normalize to the max of each img).');
end
alphaImg = photons./(maxVal*bright);
alphaImg(isnan(ccv)) = 0;
alphaImg(alphaImg > 1) = 1;

%plot the 2 images
i1 = image(photons,'Parent',axList(1),'CDataMapping','scaled');
i2 = image(ccv,'Parent',axList(2),'AlphaData',alphaImg,'CDataMapping','scaled');

%scale the images
colormap(axList(1),gray);
caxis(axList(1),[1 maxVal]);
colormap(axList(2),cmapFlip);
colorbar(axList(2),'TickLength',0);
caxis(axList(2),ccvRange);
pbaspect(axList(1),[1 size(ccv,1)/size(ccv,2) 1]);
pbaspect(axList(2),[1 size(ccv,1)/size(ccv,2) 1]);
linkprop([axList(1) axList(2)],'Position');

%now do this again without color bar to save image as tiff
f2 = figure('Position',[100 100 500 500*size(ccv,1)/size(ccv,2)]);
axList(3) = gca;
set(axList(3),'Visible','off','YLim',[0.5 0.5+size(ccv,1)],'XLim',[0.5 0.5+size(ccv,2)]);
pbaspect([1 size(ccv,1)/size(ccv,2) 1]);
hold on;
axList(4) = axes;
pbaspect([1 size(ccv,1)/size(ccv,2) 1]);
set(axList(4),'Visible','off','YLim',[0.5  0.5+size(ccv,1)],'XLim',[0.5  0.5+size(ccv,2)]);
set(axList(3),'Units','Normalized','Position', [0 0 1 1]);
set(axList(4),'Units','Normalized','Position', [0 0 1 1]);
hold on;
i3 = image(photons,'Parent',axList(3),'CDataMapping','scaled');
i4 = image(ccv,'Parent',axList(4),'AlphaData',alphaImg,'CDataMapping','scaled');
%scale the images
colormap(axList(3),gray);
caxis(axList(3),[1 maxVal]);
colormap(axList(4),cmapFlip);
caxis(axList(4),ccvRange);

fImgOnly = getframe(f2);

end

