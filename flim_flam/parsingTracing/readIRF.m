%This script reads measured IRFs in .sdt or .bin format using BioFormats.

% The script also writes a .csv file with
%the IRFs (to the directory where the .sdt files were) so that they can be
%pasted into SPC Image from the clipboard. For some reason, copying from a Matlab array to
%SPCImage directly crashes SPCImage. Column names in the .csv file
%correspond to the first filename that went into the outputted IRF
%(a relevant distinction when two IRFs are averaged in 'paired' mode).

%Dependencies: BioFormats Matlab package

%Input parameters: 
    %1. mode - either 'paired' to average IRFs in sets of 2 or 'single' to
    %use individual ones.
    %2. timeBins (format [XX XX]) - in units of ADC resolution, the time
    %bins of the input data that correspond to the IRF. For 'standard' NaI
    %fluorescein IRFs, it would be [26 36]).
    %3. outputName - base file name for saving outputs to (outputs will be
    %saved to the same folder that contain the selected .sdt files).
    %4. (optional) a path to a directory in the file system that will be a
    %convenient starting place to look for the IRFs.
    
    %sample function call: 
    %ar0726 = irfImporter('paired',[26 50],'2018-07-26_AR_irfs');
    %ar0726 = irfImporter('paired',[26 50],'2018-07-26_AR_irfs','C:\Me\SavePath');
    
%Nuts and bolts: After the IRFs are read in, they are background subtracted
%using the last quarter of the time bins, starting 50 time bins back from the end. They are then
%normalized to 15000 as the max. If mode is set to 'paired', the two IRFs
%are averaged (paired in the order that they exist in the Windows file
%system. IRFs are then cropped to the selected time bins (all other
%time bins are set to zero).

%Written by Julia Lazzari-Dean. Last edits October 12, 2019
%mode is 'single' for one IRF (no averaging) or 'paired'

function [irfOutput] = readIRF(mode,timeBins,outputName,varargin)

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
[fNames, pName] = uigetfile([pathToOpen,'*.sdt;*.bin'],'MultiSelect','on','Select raw .sdt or .bin file(s) containing IRFs.');

%figure out how many IRFs we will be processing
if(iscell(fNames))
    nIRFs = size(fNames,2);
else
    nIRFs = 1;
end

%load the first IRF to get the ADC resolution
if(nIRFs == 1)
    tempData = bfopen([pName fNames]);
else
    tempData = bfopen([pName fNames{1,1}]);
end
ADCres = size(tempData{1,1},1);

%create a struct for the raw IRFs and add the first IRF to it.
rawIRFs = zeros(ADCres,nIRFs);
for i = 1:ADCres
    if (i > size(tempData{1,1},1))
        break; %currently a hack because the .bin IRFs are coming in as varying lengths
    end
    if(size(tempData{1,1}{i,1},1) > 1) %if the IRF is an image rather than a single
        rawIRFs(i,1) = sum(sum(tempData{1,1}{i,1},1),2);
    else
        rawIRFs(i,1) = tempData{1,1}{i,1};
    end
end

%load the files in the directory
if (nIRFs > 1)
    for j = 2:nIRFs
        %open the sdt files using the Yasuda lab's parser
        try
            filePath = [pName fNames{1,j}];
            currentIRF = bfopen(filePath); % then load data
        catch % cannot load file or some other reasons an error is raised, do this
            error(['Please check the content of ' fNames{1,j}]);
        end
        
        %read the photon count data into the rawIRFs file
        for i = 1:ADCres
            if (i > size(currentIRF{1,1},1))
                break; %currently a hack because the .bin IRFs are coming in as varying lengths
            end
            if(size(currentIRF{1,1}{i,1},1) > 1) %if the IRF is an image rather than a single
                rawIRFs(i,j) = sum(sum(currentIRF{1,1}{i,1},1),2);
            else
                rawIRFs(i,j) = currentIRF{1,1}{i,1};
            end
        end
    end
end

%scale the IRFs to 15000 max (save them over into a new array)
normalizedIRFs = zeros(ADCres,nIRFs);
bkgdSubIRFs = zeros(ADCres,nIRFs);
bkgdEnd = ADCres - 50;
bkgdStart = bkgdEnd - round(ADCres/4);
for i = 1:nIRFs
    %remove the baseline counts (based on the end of the IRF)
    background = mean(rawIRFs(bkgdStart:bkgdEnd,i));
    bkgdSubIRFs(:,i) = rawIRFs(:,i) - background;
    %0 is the smallest that any value should be
    bkgdSubIRFs(bkgdSubIRFs < 0) = 0;
    
    %find the maximum and scale it to 15000
    maxVal = max(bkgdSubIRFs(:,i));
    for j = 1:ADCres
        normalizedIRFs(j,i) = bkgdSubIRFs(j,i)*15000/maxVal;
    end
end

%average (or don't) the IRFs based on the input parameters
if(strcmp(mode,'paired'))
    %average sequential IRFs in the input file for exporting (since data
    %were acquired in pairs)
    if(mod(nIRFs,2) ~= 0)
        error('Please enter an even number of IRFs for paired analyses.');
    end
    
    finalIRFs = zeros(ADCres,nIRFs/2);
    for i = 1:2:(nIRFs-1)
        for j = 1:ADCres
            %sum the two sequential values and divide by 2
            avgVal = (normalizedIRFs(j,i) + normalizedIRFs(j,i+1))/2;
            %round the nearest integer
            finalIRFs(j,(i+1)/2) = round(avgVal);
        end
    end
elseif(strcmp(mode,'single'))
    %round to the nearest integer but don't average
    finalIRFs = round(normalizedIRFs);
else
    error('IRF import mode should be either paired or single (entered as a string)');
end

%set up a struct that links file names to the data
irfOutput = struct('fName',"",'IRF',0,'directory',"");

%export the results to this struct
for ind = 1:size(finalIRFs,2)
    irfOutput(ind,1).fName = "";
    irfOutput(ind,1).directory = pName;
    if(strcmp(mode,'paired'))
        %add two files to the list
        irfOutput(ind,1).fName(1,1) = fNames{1,2*ind-1};
        irfOutput(ind,1).fName(2,1) = fNames{1,2*ind};
    elseif(strcmp(mode,'single'))
        %add one file to the list
        if(nIRFs > 1)
            irfOutput(ind,1).fName(1,1) = fNames{ind};
        else
            irfOutput(ind,1).fName(1,1) = fNames;
        end
    end
    %include the IRF only for the desired time bins
    ll = timeBins(1,1);
    ul = timeBins(1,2);
    for j = 1:ADCres
        if(j > ll && j < ul)
            irfOutput(ind,1).IRF(j,1) = finalIRFs(j,ind);
        else
            irfOutput(ind,1).IRF(j,1) = 0;
        end
    end
end

%add the word IRF to the outputname passed to the function
outputFileName = [outputName '_irfs'];

%save the struct for later use in the same folder that the IRFs are in
fullOutput = strcat(pName,outputFileName);
save(fullOutput,'irfOutput');
save(fullOutput,'mode','-append');
save(fullOutput,'timeBins','-append');
save(fullOutput,'ADCres','-append');
save(fullOutput,'rawIRFs','-append');
save(fullOutput,'normalizedIRFs','-append');
save(fullOutput,'finalIRFs','-append');

end