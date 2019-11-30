%This function is intended to be called within traceMembranes to ask the
%user to select ROIs to merge.

%The function requires an array of masks (masks), a connectivity map (cc),
%the initial image (normI), and a boolean indicating if the user has been
%trying to do this before (retry). It returns the array toMerge, which has
%a row corresponding to each group of ROIs to be merged, as well as the
%boolean mergeError indicating if the calling function should ask the user
%to try again.

%The errors that this function checks for are: less than 2 ROIs in the
%merge group or an ROI included in two different merge groups.

%Written by Julia Lazzari-Dean, October 3, 2019

function [toMerge,mergeError] = findROIsToMerge(masks,cc,normI,retry,displayScale)

if(retry)
    figTitle = 'Error - try again. Draw Around ROIs to Merge. Click to Finish.';
else
    figTitle = 'Draw Around ROIs to Merge. Click to Finish.';
end

%now ask the user to select any ROIs to merge
%convert the connectivity map to a label matrix and display it
roiMatrix2 = labelmatrix(cc);
Lrgb2 = label2rgb(roiMatrix2, 'jet', 'w', 'shuffle');
%%display the updated image on a new figure
f2 = figure('Name',figTitle,'units','normalized','outerposition',[0 0 1 1]);
%display the original image
imgH = imshow(normI,'InitialMagnification',displayScale);
hold on
%display the rois on top
himage = imshow(Lrgb2,'InitialMagnification',displayScale);
axParent = gca;
himage.AlphaData = 0.3;
text(size(normI,2)/2,0,'MERGING: Press any key to finish. To merge ROIs, click once before tracing each ROI.',...
'HorizontalAlignment','center','VerticalAlignment','bottom');

%keep finding freehand objects until the user clicks.
roiSize = 10;
numMerges = 0;
%a key press will keep the loop from being entered.
k = waitforbuttonpress;
while (k ~=1)
    %let the user create an imfreehand object
    h = imfreehand(axParent);
    %make the ROI into a mask and add it to the array
    mergeMasks(:,:,numMerges+1) = h.createMask(himage);
    %figure out how large the ROIs is and increment the number of merges.
    roiSize = sum(sum(mergeMasks(:,:,numMerges+1)));
    numMerges = numMerges+1;
    k=waitforbuttonpress;
end

%toMerge is a list of indices of ROIs to merge
toMerge = zeros(numMerges,1);
for i = 1:numMerges
    cMerge = 1;
    for j = 1:cc.NumObjects
        %iterate over all of the ROIs and see if they touch the trace
        maskTest = immultiply(masks(:,:,j),mergeMasks(:,:,i));
        if(sum(sum(maskTest)) > 0)
            %if there is any overlap between the two masks, add the index
            %to the list
            toMerge(i,cMerge) = j;
            cMerge = cMerge + 1;
        end
    end
end

%set the mergeError parameter to no errors as default, then check a few
%common error cases
mergeError = 0;
%check to see that there are no duplicates in the merge array
numROIs = nnz(toMerge);
numUniqueROIs = nnz(unique(toMerge));
if(numROIs ~= numUniqueROIs)
    %an ROI is in multiple merges
    mergeError = 1;
    disp('An individual ROI was included in more than one merge.');
end

%check that each merge has at least 1 index in it
for i = 1:numMerges
    if(nnz(toMerge(i,:)) < 2)
        mergeError = 1;
        disp('Only one ROI was identified in at least one merge.');
    end
end

end