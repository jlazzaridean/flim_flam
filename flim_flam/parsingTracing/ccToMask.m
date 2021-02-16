%This function is intended to be called within traceMembranes to process
%linear indices of an image into a binary mask that can be used to isolate
%particular regions of an image.

%It converts linear indices present in the connectivity map cc to masks
%that can be applied to an image (sampleImg is needed for sizing).

%Written by Julia Lazzari-Dean August 28, 2018

function [masks] = ccToMask(cc,sampleImg)

%%convert these ROIs into masks that can be applied to the lifetime image
%set up an object for the masks
masks = zeros(size(sampleImg,1),size(sampleImg,2),cc.NumObjects);
%fill each mask based on the values
for i=1:cc.NumObjects
    %find the linear indices for the particular ROI inside the cc struct
    idxROI = cc.PixelIdxList{1,i};
    %for each of these lists, convert to 2d indices and modify the
    %appropriate mask
    for j = 1:length(idxROI)
        [x_j, y_j] = ind2sub(size(sampleImg),idxROI(j));
        masks(x_j, y_j, i) = 1;
    end
end

end