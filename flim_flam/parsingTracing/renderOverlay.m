%Renders an overlay of a fit parameter (ccv) over the photons array
%(photons). A handle to the graphics object is returned.

%The photons array will be scaled according to photonsMax
    %photonsMax = -1 -> scale LUTs to the maximum for the photon array
    %photonsMax = any positive integer -> scale photons to this as maxium
    %for consistent display of images in a set.
    
%The dimensions of photons and ccv must match.

%The color coded value (lifetime parameter) will be rendered in the jet colormap across
%ccvRange, which should be specified as [limLow limHigh]
    
%The figure will display the title figTitle so that it is clear which image
%it corresponds to (especially important when it is called from
%renderImages or within the pxwiseAnalysis).

function [figHandle] = renderOverlay(photons,ccv,photonsMax,ccvRange,figTitle)

bright = 0.7; %enhance photon brightness a bit to make the image clearer

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
colormap(axList(2),flipud(jet));
caxis(axList(2),ccvRange);
colorbar(axList(2),'TickLength',0);
pbaspect(axList(1),[1 1 1]);
pbaspect(axList(2),[1 1 1]);
linkprop([axList(1) axList(2)],'Position');

end

