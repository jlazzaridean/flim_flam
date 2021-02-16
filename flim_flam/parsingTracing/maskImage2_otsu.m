%%This function thresholds the images using Otsu's method and removes any ROIs with
%%areas smaller than smallestObject (in units of pixels). It then uses
%%bwconncomp to assess connectivity and takes each connected group of cells
%%as one ROI.

%This version of the code implements different pre-processing depending on
%the cell type used, including sharpening the starting image various
%amounts to try to recover membranes from debris.

%Last edits 9/1/18 Julia Lazzari-Dean

function [cc] = maskImage2_otsu(normI,cellType,smallestObject)

%whether or not to run hbreak on the histogram later
spur = 0;

%based on what the image is of, pre-process the image before thresholding
if(strcmp(cellType,'A431'))
    %sharpen the image
    normI = imsharpen(normI,'amount',2,'radius',1.5);
    %adjust the smallest object threshold after sharpening - seems
    %necessary to avoid loss of signal pixels
    smallestObject = smallestObject/2;
    %best number of levels for the Otsu histogram
    nLevels = 1;
elseif(strcmp(cellType,'HEK293T'))
    %don't sharpen the image
    nLevels = 1;
elseif(strcmp(cellType,'MDA-MB-231'))
    %sharpen the image a lot to try to get thin memmbranes to stick out
    normI = imsharpen(normI,'amount',1,'radius',1.5);
    normI = imsharpen(normI,'amount',4,'radius',3);
    %break the histogram up into 2 levels
    nLevels = 2;
    %spur a little later
    spur = 1;
elseif(strcmp(cellType,'trans293T'))
    %allow more levels of relevant cells in transfected HEK293T (some dim
    %transfected cells were being removed)
    nLevels = 2;
elseif(strcmp(cellType,'MCF-7'))
    %don't sharpen the image  - HEK293T settings worked well in MCF7s
    nLevels = 1;    
else
    disp('Unrecognized cell type. Using HEK293T defaults');
    nLevels = 1;
end

%threshold the processed image
levels = multithresh(normI,nLevels);
BWA = imbinarize(normI,levels(1));

%remove small objects from the image before trying to reconnect membranes
BWA1 = bwareaopen(BWA,smallestObject);

%spur if the cell type requires it
if (spur)
    BWA1 = bwmorph(BWA1,'spur',3);
end

%process the ROIs into a connectivity matrix
cc = bwconncomp(BWA1);

end