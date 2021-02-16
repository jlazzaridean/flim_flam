%This function exports TIFF images of the total photon count present in .sdt or .bin files.

%The function prompts the user for a set of files, opens the .sdt or .bin
%files that were selected, and saves them as .TIFF into a subfolder of the
%original filepath with "_photons" added to the end of each image.

%Written by Julia Lazzari-Dean, October 11, 2019

%Dependencies: writeTIFF.m

function [] = renderPhotonsImg()

%in case the user hasn't already, add the needed subpath
addpath('parsingTracing');

%ask the user for the input files
data = readTCSPC();
pName = data(1,1).pName;
nFiles = size(data,1);

%create an output directory
oPath = strcat(pName,'rawPhotonExport');
mkdir(oPath);

for i = 1:nFiles
    oName = strcat(data(i,1).fName(1:end-4),'_photons');
    writeTIFF(sum(data(i,1).tcspc,3),oPath,oName);
end

end

