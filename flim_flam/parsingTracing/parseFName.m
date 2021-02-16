%Parses TCSPC filenames to pull out cID, imageID, repID, and frameID.

%cID - identifier for sample or coverslip
%imageID - identifier for each unique field of view
%repID - identifier for the number of times this FOV was recorded from
%frameID - identifier for the frame index in a time series 

%The function parses the file name strings to determine a coverslip and
%image number for each file. It assumes that the file names are in the form
%   YYYY-MM-DD_CC-II-RR_otherStuff_XXX 
%All numerical fields can be varying lengths as long as they are separated
%by the expected delimiter. The underscore before the frameID is not
%required as long as the preceding character is not a number.

%If a repID RR is not assigned, it will default to 1 (the first image on
%this field of view).

%Mode is whether the image is the only image in the series for that
%cID_imageID or if it is part of a time series. expects values of 'individ'
%or 'series'. If the mode is 'series', the function will look at the end of
%the filename and return the number that is there (will find a variety of
%lengths). If the mode is 'individ', the frameID returned is -1 and the end
%of the string is not parsed.

%Last edited October 10, 2019 by Julia Lazzari-Dean
%Dependencies: parseIDString

function [date, cID, imageID, repID, frameID] = parseFName(fName,mode)

%take off the .sdt file extension
baseName = fName(1:end-4);
date = baseName(1:10);
frameID = -1;

%start parsing the base file name to identify the coverslip and image IDs
%remove the date from the name to get to the coverslip ID
string1 = baseName(12:end);

%remove the word 'slip' if it was included
if(strcmp(string1(1:4),'slip'))    
    string1 = string1(5:end);
end

[cID,imageID,repID] = parseIDString(string1);

%if this image is part of a series, find a frame ID
if(strcmp(mode,'series'))
    test = 1;
    readCount = 0;
    nEntries = size(string1,2);
    while (~isnan(test))
        if readCount > nEntries
            error('Missing delimiter; cannot isolate frame ID');
        end
        entry = string1(end-readCount);
        test = str2double(entry);
        readCount = readCount + 1;
    end
    frameID = str2double(string1((end+2-readCount):end));
end

end

