%%This script produces membrane ROIs based on a photon count image
%%"photons" and returns these ROIs as binary masks. The goal is to identify
%%free-standing cell groups for analysis of FLIM-based membrane potential.

%%It uses the otsu version of create mask to find objects in the image.
%%This function calls v2 of that script, which has cell type specific
%%image processing implemented for A431, CHO, HEK293T, MCF-7, and MDA-MB-231
%%cells.

%%This function supports debris removal and merging of ROIs that were split
%%up but shouldn't have been. The script first asks the user to remove any
%%debris (by right clicking) and then invites the user to merge any ROIs by
%%circling the ROIs to merge. It then returns the lifetime results in each
%%ROI, as well as a boolean as to whether the user was satisfied with the
%%final ROIs.

%%Last edited October 3, 2019 by Julia Lazzari-Dean

function [masksFinal,goodROIs] = traceMembranes4(photons, cellType, smallestObject)

%figure out the magnification to display the image
nativeWidth = size(photons,2);
set(0,'units','pixels');
Pix_SS = get(0,'screensize');
width_SS = Pix_SS(1,3);
scale = width_SS/nativeWidth * 0.4 * 100;

%normalize the photons image for display and masking
%divide everything by the largest value in the image
large = max(photons);
largest = max(large);
%normalize the image to this largest value
normI = photons./largest;

%mask the image using the desired method & return a connectivity matrix
cc = maskImage2_otsu(normI,cellType,smallestObject);

%set up the rois for easy viewing by the user
%convert the connectivity map into a label matrix
roiMatrix = labelmatrix(cc);
%color code the label matrix
Lrgb = label2rgb(roiMatrix, 'jet', 'w', 'shuffle');

%%Prompt the user to remove any large debris
%open a full-screen figure the user can interact with
f1 = figure('Name','Left Click Regions to Remove. Enter to Finish.',...
    'units','normalized','outerposition',[0 0 1 1]);
%display the original image
imshow(normI,'InitialMagnification',scale);
hold on
%display the rois on top
himage = imshow(Lrgb,'InitialMagnification',scale);
himage.AlphaData = 0.3;
text(size(normI,2)/2,0,'CLEANING: Left click ROIs (debris) to remove. Enter when finished.',...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

%record the areas the user got rid of
[x,y] = getpts(f1);
%round to full pixels
x = round(x);
y = round(y);
%find the index of the pixel that was selected
linInd = sub2ind(size(normI),y,x);

%search to find the regions that contain these indexed points
remove = zeros(1,length(linInd));
for i = 1:length(linInd)
    for j=1:length(cc.PixelIdxList)
        if (ismember(linInd(i), cc.PixelIdxList{j}))
            remove(i) = j;
        end
    end
end
%return some sort of error if the region hasn't been assigned yet
if(ismember(0,remove))
   disp('Some identified pixels were not in an ROI. Ignoring these.'); 
   where_zero = find(remove == 0);   
   remove(where_zero) = [];
end

%remove duplicates in case the same ROI was clicked twice
remove = unique(remove);

%sort the regions to be removed into ascending order so that the 'remove
%from a list as you iterate over it' bit works
remove = sort(remove);

%remove regions that contain these points, adjusting indices as you go
numRemoved = 0;
for i = 1:length(remove)
    if((remove(i)-numRemoved) < 1)
        disp('What is going on');
    end
    cc.PixelIdxList(remove(i)-numRemoved) = [];
    cc.NumObjects = cc.NumObjects - 1;
    numRemoved = numRemoved+1;
end

%convert these ROIs into masks that can be applied to the lifetime image
masks = ccToMask(cc,photons);

mergeError = 1;
nAttempts = 0;
while(mergeError == 1)
    %check if this loop has run before so that an error message can be run
    %for the user
    if (nAttempts == 0)
        retry = 0;
    else
        retry = 1;
    end
    %to merge is a matrix of mask indices to merge
    %each row is a separate set of indices to merge
    %columns correspond to the mask index
    [toMerge,mergeError] = findROIsToMerge(masks,cc,normI,retry,scale);
    nAttempts = nAttempts + 1;
end

%actually merge the ROIs that were returned in the original connectivity
%matrix
%set up logical indexing to remove the merged indices
saveIndex = true(1,cc.NumObjects);
for i = 1:size(toMerge,1)
    %iterate over all of the rows of toMerge and merge the ROIs
    %save the first ROI in each list - so j starts at 2
    for j = 2:size(toMerge,2)
        if(toMerge(i,j) ~= 0)
            %for each group, add all of the ROIs to the first one
            origList = cc.PixelIdxList{1,toMerge(i,1)};
            addendum = cc.PixelIdxList{1,toMerge(i,j)};
            cc.PixelIdxList{1,toMerge(i,1)} = cat(1,origList,addendum);
            %set the logical index to false so it is deleted later
            saveIndex(1,toMerge(i,j)) = 0;
        end
    end    
end

%remove the indices that are now duplicate
ccMerged = cc;
for i=size(ccMerged.PixelIdxList,2):-1:1
    %iterate backwards to avoid issues with indices changing
    if(saveIndex(1,i) == 0)
        ccMerged.PixelIdxList(:,i) = [];
    end
end
%update the number of objects in the label matrix
ccMerged.NumObjects = sum(saveIndex);

%close the other figures to get ready to display the final ROIs
close all;
%convert the updated list to a new color coded label matrix
roiMatrix3 = labelmatrix(ccMerged);
Lrgb3 = label2rgb(roiMatrix3, 'jet', 'w', 'shuffle');
%%display the updated image on a new figure
f3 = figure('Name','Final ROIs','units','normalized','outerposition',[0 0 1 1]);
%display the original image
imshow(normI,'InitialMagnification',scale);
hold on
%display the rois on top
himage = imshow(Lrgb3,'InitialMagnification',scale);
himage.AlphaData = 0.3;
text(size(normI,2)/2,0,'FINALIZE: Press Any Key to Approve. Click to Reject & Re-Do.',...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

%set the goodROIs boolean based on the user's approval (or lack thereof)
goodROIs = waitforbuttonpress;

%convert these ROIs into masks that can be applied to the lifetime image
%set up an object for the masks
masksFinal = ccToMask(ccMerged,photons);

close all;

end
