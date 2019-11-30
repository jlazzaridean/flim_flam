%Reads a .csv of metadata linking ROIs to images and generates a structure.

function [roiMD, roiMDName] = readROIMD(varargin)

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

%prompt the user for a config file .csv (basically TCSPC system parameters).
%read this information into a structure that can be parsed by the other
%scripts so that all variables don't have to be entered individually
[roiMDName,pName] = uigetfile([pathToOpen '*.csv'],'Please select .csv ROI METADATA file');
fullPath = [pName fsep roiMDName];

T = readtable(fullPath);
roiMD = table2struct(T);

end

