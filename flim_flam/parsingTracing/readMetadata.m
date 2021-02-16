%Reads a .csv containing metadata about the recordings and generates a structure. 

%It also returns the file name of the metadata so that the user can approve
%it.

function [md,fName] = readMetadata(varargin)

%open the convenient path if it is provided to the function
fsep = filesep;
if(size(varargin,2) == 1)
    pathToOpen = varargin{1,1};
else
    pathToOpen = pwd;
end
if(pathToOpen(end) ~= fsep)
    pathToOpen(end+1) = fsep;
end

[fName,pName] = uigetfile([pathToOpen '*.csv'],'Please select .csv METADATA file');
fullPath = [pName fsep fName];

T = readtable(fullPath);
md = table2struct(T);
fields = fieldnames(md);

%make sure cIDs and IRFindex were specified
if(~any(strcmp(fields,'irf_index')))
    error('Please specify the parameter IRFindex in the metadata');
elseif(~any(strcmp(fields,'cID')))
    error('Please specify the parameter cID in the metadata');
end


end

