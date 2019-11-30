%Reads a user-specified set of .TIFF files in and interprets them as binary masks.

%An optional second input parameter specifies the directory to open first
%when asking the user to specify ROI filenames.

%This function operates in two modes (specified by mode); more details on
%each are available in the manual.
%   'csv' - Ask for a .csv input file that specifies ROI names and the
%   various identifiers (cID, imageID, repID, frameID...as relevant). These
%   ROI names are then used to match ROIs with files. There must be at
%   least one metadata entry *per file.* 
%   'parse' - Use the structure of the ROI names to infer which image files
%   they should be applied to. For this to work, the ROIs must be named
%   regularly, where CC = cID, II = imageID, RR = repID, ZZ = ROI number
%   for that image: CC-II-RR_ZZ (missing RR will be ignored; underscore
%   separates ROI number from other identifiers)
%   PARSE DOES NOT SUPPORT USING DIFFERENT ROIS ON
%   DIFFERENT FRAMES OF THE SAME TIME SERIES. The same parsed ROIs will be
%   applied to every frame of the time series.

%Written by Julia Lazzari-Dean, last edits October 12, 2019.

%Dependencies: readROIMD, parseIDString

function [rois] = readImageJROIs(mode,varargin)

%get current working directory and system file separator
fsep = filesep;
%open the convenient path if it is provided to the function
if(size(varargin,2) == 1)
    pathToOpen = varargin{1,1};
else
    pathToOpen = pwd;
end
if(pathToOpen(end) ~= fsep)
    pathToOpen(end+1) = fsep;
end

%ask for the rois
[fNames, pName] = uigetfile([pathToOpen,'*.tif;*.tiff'],'MultiSelect','on','Select exported TIFF ROIs.');

if strcmp(mode,'csv')    
    %find the ROI file extension for later matching
    if size(fNames,2) == 1
        fNameSample = fNames;
    else
        fNameSample = fNames{1,1};
    end
    extInd = strfind(fNameSample,'.');
    fileExt = fNameSample(extInd:end);
    
    %ask for and get the csv metadata if applicable
    [rois,roiMDName] = readROIMD(pathToOpen);
    mdFields = fieldnames(rois);
    defaultRepID = 0;
    defaultFrameID = 0;
    if (~any(contains(mdFields,'cID')) || (~any(contains(mdFields,'imageID'))))
        error('Provided ROI MD is missing cID and/or imageID information');
    end
    if ~any(contains(mdFields,'repID'))
        %if replicate ID is not provided, fill with 1 later on
        defaultRepID = 1;
    end
    if ~any(contains(mdFields,'repID'))
        %if frame ID is not provided, fill with 1 later on
        defaultFrameID = 1;
    end
    nROIs = size(rois,1);
    nFiles = size(fNames,2);
    
    %open all of the ROIs and read in their masks
    roisFromFiles(1:nFiles,1) = struct('roiName',0,'masks',0);
    for i=1:nFiles
        if nFiles == 1
            roisFromFiles(1,1).roiName = fNames;
        else
            roisFromFiles(i,1).roiName = fNames{1,i};
        end        

        fullROIPath = strcat(pName,fsep,roisFromFiles(i,1).roiName);
        t = Tiff(fullROIPath,'r');
        img = read(t);
        img(img>0) = 1; %make binary
        roisFromFiles(i,1).masks = img;
    end
        
    %iterate over all of the ROIs and find the corresponding masks
    %there may be more than one metadata entry for a given ROI if it is
    %supposed to be applied to various images.    
    for i=1:nROIs
        rois(i,1).roiPath = pName;
        rois(i,1).roiMDName = roiMDName;        
        
        %add a file extension and look for matches
        testName = strcat(rois(i,1).roiName,fileExt);
        for j=1:nFiles
            if strcmp(testName, roisFromFiles(j,1).roiName)
                myROI = roisFromFiles(j,1);
                break;
            end
        end
        
        %add the matched data into the rois structure
        if isequal(size(myROI),[1 1]) %unique match was found
            rois(i,1).masks = myROI(1,1).masks;
        else
            error('A unique match for at least one roiName in metadata was not found in provided ROIs.');
        end
        
        if defaultFrameID %fill in the defaults on frameID and repID if they were not in metadata
            rois(i,1).frameID = 1;
        end
        if defaultRepID
            rois(i,1).repID = 1;
        end
    end
    
elseif strcmp(mode,'parse')
    nROIs = size(fNames,2);
    rois(1:nROIs,1) = struct('roiPath',0,'roiName',0,'masks',0,'cID',0);
        
    for i=1:nROIs
        %save the file name and path name of the ROIs
        rois(i,1).roiPath = pName;
        if nROIs == 1
            rois(1,1).roiName = fNames;
        else
            rois(i,1).roiName = fNames{1,i};
        end
        
        %find the name of this ROI
        extInd = strfind(rois(i,1).roiName,'.');
        name = rois(i,1).roiName(1:extInd-1);
        if contains(name,'_')
            %if there are multiple ROIs for cID/imageID/repID/frameID set
            splitInd = strfind(name,'_');
            string1 = name(1:splitInd-1);
            roiIDString = name(splitInd+1:end);
            rois(i,1).roiID = str2double(roiIDString);
            [cID,imageID,repID] = parseIDString(string1);
        else
            %if there is only one ROI for cID/imageID/repID/frameID set
            rois(i,1).roiID = 1;
            [cID,imageID,repID] = parseIDString(name);
        end
        %assign other parameters
        rois(i,1).cID = cID;
        rois(i,1).imageID = imageID;
        rois(i,1).repID = repID;
        
        %if anything didn't assign correctly, throw an error
        if (any(isnan(rois(i,1).cID)) || any(isnan(rois(i,1).imageID)) || any(isnan(rois(i,1).repID)))
            error(['Incorrect parsing of roi filename ' rois(i,1).roiName]);
        end
        
        %now actually read the tiff file and save image matrix as mask
        fullROIPath = strcat(pName,fsep,rois(i,1).roiName);
        t = Tiff(fullROIPath,'r');
        img = read(t);
        img(img>0) = 1; %make binary
        rois(i,1).masks = img;
    end
else
    error('ROI mode must be either parse or csv.');
end

end

