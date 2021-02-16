%This is the main function to import .bin or .sdt data and generate decays to fit. 

%After recording information from the user (see
%below), it will return three structures and save versions of these three
%structures to .csv and .mat in the folder where the data were contained. 
    %1. A structure of data, paired with metadata, with decays to be fit.
    %2. The unpaired metadata structure for review.
    %3. A configuration structure, which is a required input for all of the
    %exponential fitting functions.

%Written by Julia Lazzari-Dean, October 2019. No guarantees or support are
%provided or implied. More information is provided in the manual.

%Sample function call:
    %[myData,myMetadata,myConfig] = importData('single','sampleOutputName');

%The function will prompt the user for at least four different things in 
%the process of importing data. More information about requirements for
%each of these things is included in the user manual.
    %1. Config file (*.csv) with the system parameters and fit parameters. Samples
    %are included with this code.
    %2. Photon histogram files (.bin or .sdt) containing the data to be
    %analyzed.
    %3. Metadata file (*.csv) with information about each record. The
    %metadata must contain the field cID (coverslip/sample identifier) that
    %matches identifiers found in the input data (see below), as well as
    %the field irf_index, which specifies which IRF number will be used for
    %processing the files. It may contain imageID, replicateID, or frameID
    %information. Other fields will be carried forward with unique
    %cID, imageID, replicateID, frameID sets.
    %4. Measured IRF files (.bin or .sdt) containing IRF for the data,
    %either as a single decay 1 x N or as an image X by Y by N. IRFs will
    %be reshaped a single N x 1 decay per file or file pair (if 'paired' is
    %specified and IRF averaging is desired).

%The way the photon histograms are converted to decays depends on the
%parameter importMode. The following four options are acceptable values of
%importMode:
    %%%mode 'single' takes images or single decays and collapses it down so
    %that one decay is generated per image    
    %%%mode 'pxwise' keeps the decays as an X by Y by N array, where X by Y
    %is the spatial resolution of the input image and N is the ADC
    %resolution
        %you will have an opportunity to bin later during fitting but not
        %here    
    %%%mode 'global_btm' is the start of a global analysis (one fit per ROI)
    %using ROIs defined by thresholding in the "batch trace membranes" code
        %to use this mode, the metadata file must contain the following two
        %fields (exact match): cellType, smallestObject    
    %%%mode 'global_fiji' is the start of global analysis (one fit per ROI) that
    %will prompt the user for a .zip file containing FIJI ROIs to apply to
    %the image
        %There are two 'sub-modes' here, which the software will ask the
        %user to dictate. The first option is to parse ROI file names to extract
        %information about the correct image match. The second option is to
        %provide a separate ROI matching metadata file. Additional details
        %about formatting and requirements for both modes are available in
        %the user manual.
              
%Dependencies (all .m): ccToMask, findROIsToMerge, makeSingleDecay,
%maskImage2_otsu, maskRawData, parseFName, parseIDString, readBin,
%readConfig, readImageJROIs, readIRF, readMetadata, readROIMD, readTCSPC,
%traceMembranes4, writeTIFF

function [dataToSave,md,config] = importData(importMode,oName)

%add subfolders to the path
addpath('parsingTracing');

%one of the only hardcoded things - the convenient path on my computer
convenientPath = 'C:\Users\JLD\PhD\ImagingData\';
fsep = filesep;

%get the configuration file
[config,configName] = readConfig(convenientPath);

%ask the user for the input files
data = readTCSPC(convenientPath);
pName = data(1,1).pName;
nFiles = size(data,1);

%ask the user for the metadata (default to same path as images)
[md, mdName] = readMetadata(pName);
mdFields = fieldnames(md);
readMetadataRepID = 0;
readMetadataImageID = 0;
readMetadataFrameID = 0;
isTimeSeries = 0;
if(any(strcmp(mdFields,'repID')))
    readMetadataRepID = 1;
end
if(any(strcmp(mdFields,'imageID')))
    readMetadataImageID = 1;
end
if(any(strcmp(mdFields,'frameID')))
    readMetadataFrameID = 1;
end
if(any(strcmp(mdFields,'nFrames')) || any(strcmp(mdFields,'frameID')))
    isTimeSeries = 1;
end

%verify that the metadata contained the necessary parameters if mode is
%'global_btm'
if(strcmp(importMode,'global_btm'))
    if(~any(strcmp(mdFields,'smallestObject')))
        error('If using global_btm as mode, you must specify smallest_object in the metadata.');
    elseif(~any(strcmp(mdFields,'cellType')))
        error('If using global_btm as mode, you must specify cell_type in the metadata.');
    end
    %create a folder in the directory for the raw BTM masks
    maskPath = strcat(data(1,1).pName,fsep,'GlobalBTM_Masks');
    mkdir(maskPath);
end

%if we are going to process with spatial resolution, prepare to render the
%photon images
if ~strcmp(importMode,'single')
    photonPath = strcat(data(1,1).pName,'rawPhotonExport');
    mkdir(photonPath);
end

%ask the user for the IRF files
irfStruct = readIRF(config(1,1).IRFmode,config(1,1).IRFlim,oName,pName);
nIRFs = size(irfStruct,1);

%test that all of the IRFs are specified in the metadata
md_irf_list = [md(:,1).irf_index]';
irf_use_test = ismembertol(1:nIRFs,md_irf_list);
%test if there are any zeros in the array
irf_use_count = nnz(~irf_use_test);
if (irf_use_count > 0)
    disp('Warning - number of IRFs imported is more than the number of IRFs referred to in the metadata.');
end

%find the expected number of ADC channels and verify that all of the IRFs
%are the expected length (decays are verified below)
ADCres = config(1,1).ADCres;
for i=size(irfStruct,1)
    if(size(irfStruct(i,1).IRF,1) ~= ADCres)
        error('Some of the selected IRFs contain a number of ADC channels inconsistent with the config specifications.');
    end
end

%Display a bit of information about selections so that the user can see if
%they want to continue
disp("You have selected the following options for fitting in mode " + string(importMode));
disp(sprintf('\tConfig: ') + " " + string(configName));
disp(sprintf('\tData files (%d total) beginning with: ',nFiles) + " " + string(data(1,1).fName));
disp(sprintf('\tMetadata: ') + " " + string(mdName));
disp(sprintf('\tIRF files (%d total) beginning with: ',size(irfStruct,1)) + " " + string(irfStruct(1,1).fName(1,1)));
%receive user input
answer = "";
while(~strcmp(answer,'Y') && ~strcmp(answer,'y') && ~strcmp(answer,'N') && ~strcmp(answer,'n'))
    prompt = 'Accept (Y) or decline {N): ';
    answer = input(prompt,'s');
end
%process the user input and either exit or continue.
if(strcmp(answer,'N')  || strcmp(answer,'n'))
    disp('File selections rejected. Please restart script to continue.');
    dataToSave = -1;
    config = -1;
    md = -1;
    return
end

%if mode is global_fiji, ask the user if they would like to parse ROI names from the
%save file or upload an roi metadata
if(strcmp(importMode,'global_fiji'))
    fijiMode = -1;
    fprintf('\n');
    disp('How would you like to match ROIs with data files?');
    while(~(fijiMode == 1) && ~(fijiMode ==2))
        prompt = 'Press 1 to parse IRF filenames or 2 to upload an ROI metadata (*.csv): ';
        fijiMode = input(prompt);
    end
    if(fijiMode == 1)
        roiS = readImageJROIs('parse',convenientPath);
    elseif (fijiMode == 2)
        roiS = readImageJROIs('csv',convenientPath);
    end
end

%match the files with metadata and create decays based on the user's
%preferences
for i=1:nFiles
    %check that all exports contain the same number of ADC channels
    if(size(data(i,1).tcspc,3) ~= ADCres)
        error('Some of the included files contain numbers of ADC channels different from the config specification.');
    end
    
    %parse the string to find cID, imageID, and frameID (if applicable)
    if isTimeSeries
        [data(i,1).date, data(i,1).cID, data(i,1).imageID, data(i,1).repID, data(i,1).frameID] = parseFName(data(i,1).fName,'series');
    else
        [data(i,1).date, data(i,1).cID, data(i,1).imageID, data(i,1).repID, ~] = parseFName(data(i,1).fName,'individ');
    end
    
    %add the correct metadata to this structure. This will always match
    %cID. If imageID and/or frameID are present in the metadata structure,
    %it will match those too.
    cIDMatch = md([md.cID]' == data(i,1).cID);

    if readMetadataImageID && readMetadataRepID && readMetadataFrameID
        %cID, image ID, replicate ID, and frame ID were all provided.
        iIDMatch = cIDMatch([cIDMatch.imageID]' == data(i,1).imageID);
        rIDMatch = iIDMatch([iIDMatch.repID]' == data(i,1).repID);
        fIDMatch = rIDMatch([rIDMatch.frameID]' == data(i,1).frameID);
        if (isempty(fIDMatch) || size(fIDMatch,1) > 1) %if some aspect didn't assign correctly
            error(strcat('The coverslip ID ',num2str(data(i,1).cID), ' and image ID ',...
                num2str(data(i,1).imageID),' did not uniquely match metadata.'));
        end
        reduced = rmfield(fIDMatch,{'cID','imageID'});
    elseif readMetadataImageID && readMetadataFrameID
         %cID, image ID, and frame ID were provided in metadata. no repID.  
        iIDMatch = cIDMatch([cIDMatch.imageID]' == data(i,1).imageID);
        fIDMatch = iIDMatch([iIDMatch.frameID]' == data(i,1).frameID);
        if (isempty(fIDMatch) || size(fIDMatch,1) > 1) %if some aspect didn't assign correctly
            error(strcat('The coverslip ID ',num2str(data(i,1).cID), ' and image ID ',...
                num2str(data(i,1).imageID),' did not uniquely match metadata.'));
        end
        reduced = rmfield(fIDMatch,{'cID','imageID'});        
    elseif readMetadataImageID && readMetadataRepID
        %cID, image ID, and replicate ID were provided in metadata. no
        %frame ID.
        iIDMatch = cIDMatch([cIDMatch.imageID]' == data(i,1).imageID);
        rIDMatch = iIDMatch([iIDMatch.repID]' == data(i,1).repID);
        if (isempty(rIDMatch) || size(rIDMatch,1) > 1) %if some aspect didn't assign correctly
            error(strcat('The coverslip ID ',num2str(data(i,1).cID), ', image ID ',num2str(data(i,1).imageID),...
                ', and rep ID ',num2str(data(i,1).repID),' did not uniquely match metadata.'));
        end
        reduced = rmfield(rIDMatch,{'cID','imageID','repID'});
    elseif readMetadataRepID
        %cID and replicate ID were provided in the metadata
        rIDMatch = cIDMatch([cIDMatch.imageID]' == data(i,1).repID);
        if (isempty(rIDMatch) || size(rIDMatch,1) > 1) %if some aspect didn't assign correctly
            error(strcat('The coverslip ID ',num2str(data(i,1).cID), ...
                ' and rep ID ',num2str(data(i,1).repID),' did not uniquely match metadata.'));
        end
        reduced = rmfield(rIDMatch,{'cID','repID'});        
    elseif readMetadataImageID
        %cID and imageID were provided in the metadata
        iIDMatch = cIDMatch([cIDMatch.imageID]' == data(i,1).imageID(1,1));
        if (isempty(iIDMatch) || size(iIDMatch,1) > 1) %if some aspect didn't assign correctly
            error(strcat('The coverslip ID ',num2str(data(i,1).cID), ' and image ID ',...
                num2str(data(i,1).imageID),' did not uniquely match metadata.'));
        end
        reduced = rmfield(iIDMatch,{'cID','imageID'});
    else
        %only the cID is provided in the metadata
        if (isempty(cIDMatch) || size(cIDMatch,1) > 1) %if the cID didn't assign correctly
            error(['The coverslip ID ' num2str(data(i,1).cID) ' did not uniquely match metadata.']);
        end
        reduced = rmfield(cIDMatch,'cID');
    end
    %take the reduced metadata structure and copy it into the data
    red_fields = fieldnames(reduced);
    for j=1:size(red_fields,1)
        data(i,1).(red_fields{j,1}) = reduced.(red_fields{j,1});
    end
    
    %if the frameID is out of the appropriate range for the metadata, quit
    if isTimeSeries && any(strcmp(mdFields,'nFrames'))
        if(data(i,1).frameID > data(i,1).nFrames || data(i,1).frameID < 0)
            error(strcat('Frame ID out of range on file entry ',num2str(i)));
        end
    end
    
    %assign the correct IRF based on the metadata
    irf_index = data(i,1).irf_index;    
    %return an error if the category was not found
    if(irf_index < 1 || irf_index > nIRFs)
        error('Could not assign an IRF based on metadata.');
    end
    data(i,1).IRF = irfStruct(irf_index,1).IRF;
    data(i,1).irfName = irfStruct(irf_index,1).fName;
    
    %now generate the decays to fit from the tcpsc array
    if(strcmp(importMode,'single'))        
        offsetProvided = any(strcmp(mdFields,'offset')); %if offset was provided as a metadata parameter, subtract it
        if offsetProvided
            offset = data(i,1).offset;
        else
            offset = 0;
        end
        data(i,1).decays = makeSingleDecay(data(i,1).tcspc,offset);        
    elseif(strcmp(importMode,'pxwise'))
        %render the photons image for each image but otherwise don't do
        %mmuch else
        writeTIFF(sum(data(i,1).tcspc,3),photonPath,strcat(data(i,1).fName(1:end-4),'_rawPhotons'));
    elseif(contains(importMode,'global'))
        %generate photons with spatial resolution and call the ROI
        %identification script
        photons = sum(data(i,1).tcspc,3);
        writeTIFF(photons,photonPath,strcat(data(i,1).fName(1:end-4),'_rawPhotons'));
        if(contains(importMode,'btm'))
            goodROIs = 0;
            while(~goodROIs)
                [masks, goodROIs] = traceMembranes4(photons,data(i,1).cellType,data(i,1).smallestObject);
            end
            data(i,1).masks = masks;
            nROIs = size(data(i,1).masks,3);
            data(i,1).roiID = zeros(nROIs,1);
            for j = 1:nROIs
                writeTIFF(data(i,1).masks(:,:,j),maskPath,strcat(data(i,1).fName(1:end-4),'_mask',num2str(j)));
                data(i,1).roiID(j,1) = j;
            end
        elseif(contains(importMode,'fiji'))
            writeTIFF(sum(data(i,1).tcspc,3),photonPath,strcat(data(i,1).fName(1:end-4),'_rawPhotons'));
            %match up the ROIs with the files as appropriate - need to
            %match coverslipID, imageID (if provided), and replicate ID
            %(should set to 1 by default so no issue if not provided)
            myROI = roiS([roiS.cID]' == data(i,1).cID);
            myROI = myROI([myROI.repID]' == data(i,1).repID);
            if (data(i,1).imageID ~= -1)
                myROI = myROI([myROI.imageID]' == data(i,1).imageID);
            end
            %if the ROIs were parsed it will by default apply to the first
            %frame (or all frames) of an images set. frame ID must also
            %match for the csv imports. frameID will default to 1 if
            %missing in ROI metadata
            if (isTimeSeries && answer == 2)
                myROI = myROI([myROI.frameID]' == data(i,1).frameID);
            end
            %copy over the matched ROI information into the data structure
            nROIs = size(myROI,1);
            if nROIs == 0
                disp("No ROIs found for record " + num2str(i));
                continue %go to next iteration; don't try to make decays
            end
            data(i,1).masks = zeros(size(data(i,1).tcspc(:,:,1)));
            data(i,1).roiName = strings(nROIs,1);
            data(i,1).roiPath = myROI(1,1).roiPath;
            data(i,1).roiID = zeros(nROIs,1);
            for j = 1:nROIs
                if isequal(size(myROI(j,1).masks),size(data(i,1).tcspc(:,:,1)))
                    data(i,1).roiName(j,1) = myROI(j,1).roiName;
                    data(i,1).masks(:,:,j) = myROI(j,1).masks;
                    data(i,1).roiID(j,1) = myROI(j,1).roiID;
                else
                    error('Dimension mismatch between TCSPC data and ROIs.');
                end 
            end
        else
            error('ROI mode for global must be either FIJI or BTM');
        end
        
        %apply masks to the data
        data(i,1).decays = maskRawData(data(i,1).tcspc,data(i,1).masks);        
    else
        error('Unrecognized import mode');
    end
end

%test that all metadata were used & warn if it wasn't
cIDs_md = [md(:,1).cID]';
cIDs_files = [data(:,1).cID]';
mdTest = ismembertol(cIDs_md,cIDs_files);
%test if there are any zeros in the array
mdTestCount = nnz(~mdTest);
if (mdTestCount > 0)
    disp('Warning - some metadata entries are unused.');
end

%save the data array to the output name in "path"
fullOName = [pName oName '_raw'];
switch importMode
    case 'single'
        dataToSave = rmfield(data,'tcspc');
        save(fullOName,'dataToSave','md','config','importMode');
        writetable(struct2table(dataToSave),[fullOName '.csv']);
    case 'global_btm'
        dataToSave = rmfield(data,'tcspc');
        save(fullOName,'dataToSave','md','config','importMode');
        writetable(struct2table(dataToSave),[fullOName '.csv']);
    case 'global_fiji'
        dataToSave = rmfield(data,'tcspc');
        save(fullOName,'dataToSave','md','config','roiS','importMode');        
        writetable(struct2table(dataToSave),[fullOName '.csv']);
    case 'pxwise'
        dataToSave = data;
        save(fullOName,'data','md','config','importMode','-v7.3');
    otherwise
        error('Unrecognized import mode');
end


end

