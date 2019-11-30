%Reads in TCSPC photon histograms in either .sdt or .bin format.

%Prompts the user to specify a list of files to read.

%Returns a structure containing the photon histogram in the field tcspc,
%along with the filename and path corresponding to each entry.

%Dependencies: BioFormats MATLAB library, readBin.m

function [rawDataStruct] = readTCSPC(varargin)

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

%ask the user to select files from a directory & verify file formats
[fNames, pName, ~] = uigetfile([pathToOpen '*.sdt;*.bin'],'MultiSelect','on','Select .sdt or .bin image file(s).');
switch class(fNames)
    case 'double'
        disp('No input file selected. Exiting.')
        return
    case 'cell'
        % Do nothing
    case 'char'
        fNames = {fNames}; % format to cell to work with stuffs
    otherwise
        error('Error.')
end

%processing names and setting up data structure for image data
nFiles = length(fNames);
rawDataStruct = struct('fName',"",'tcspc',0); %Empty struct, name holder, image data holder

%iterate over the  .sdt files and read the data into the structure
for i=1:nFiles
    %save the file name and path name to the structure
    filePath = [pName fsep fNames{1,i}];
    rawDataStruct(i,1).fName = fNames{1,i};
    rawDataStruct(i,1).pName = pName;
    
    %open the correct file parser and read data
    if (strcmp(fNames{1,i}(end-2:end),'bin'))
        rawDataStruct(i,1).tcspc = readBin(filePath);     
    else
        tempData = bfopen(filePath);
        firstImg = tempData{1,1}{1,1};
        imgStack = zeros(size(firstImg,1),size(firstImg,2),size(tempData{1,1},1));
        for k = 1:size(imgStack,3) %dimension 3 is time in ns
            imgStack(:,:,k) = tempData{1,1}{k,1};
        end
        rawDataStruct(i,1).tcspc = imgStack;
    end
end

end