%Helper function to extract numbers of varying lengths from a string of the
%format XX-XX-XX, where each XX can be any length as long as they are
%separated by dashes.

function [cID,imageID,repID] = parseIDString(string1)

%instantiate nonsensical values to help with parsing the string
imageID = -1;
repID = 1; %if no repID is found, will assume that this is the first (and only) replicate of a given FOV
imageIDExist = 1;
repIDExist = 1;

%read back from the beginning of the string to find the first dash or
%underscore
test = 1;
readCount1 = 1; %reading forward from the front of the string
nEntries1 = size(string1,2);
while (~isnan(test))
    if readCount1 > nEntries1
        readCount1 = readCount1 + 1; %increment the counter and exist as if you found a dliminter
        break;
    end    
    entry = string1(readCount1);
    test = str2double(entry);
    readCount1 = readCount1 + 1;
end
if(~strcmp(entry,'-')) %if there is anything other than a dash after the final entry
    imageIDExist = 0;
end
if(readCount1<2)
    error('No coverslip ID found after the date');
else
    cID = str2double(string1(1:(readCount1-2)));
    string2 = string1(readCount1:end);
end

%now start the search for the image ID
if ~ imageIDExist
    disp(['coverslip ID without imageID found for ' string1]);
else
    %read back from the beginning of the string to find the first dash or
    %underscore
    test = 1;
    readCount2 = 1; %reading forward from the front of the string
    nEntries2 = size(string2,2);
    while (~isnan(test))
        if readCount2 > nEntries2
            readCount2 = readCount2+1; %increment the counter as if you found a delimiter and exit loop
            break;
        end        
        entry = string2(readCount2);
        test = str2double(entry);
        readCount2 = readCount2 + 1;
    end
    if(readCount2 < 2)
        disp(['ImageID not parsed correctly for ' string1]);
    else
        imageID = str2double(string2(1:(readCount2-2)));
        string3 = string2(readCount2:end);
    end
    if(~strcmp(entry,'-')) %there is anything other than a dash next
        repIDExist = 0;
    end
    
    if(repIDExist) %if there was a dash after the image ID
        test = 1;
        readCount3 = 1; %reading forward from the front of the string
        nEntries3 = size(string3,2);
        while (~isnan(test))
            if readCount3 > nEntries3
                readCount3 = readCount3+1; %increment the counter as if you found a delimiter and exit loop
               break;
            end            
            entry = string3(readCount3);
            test = str2double(entry);
            readCount3 = readCount3 + 1;
        end
        if(readCount3 < 2)
            disp(['RepID not parsed correctly for ' string1]);
        else
            repID = str2double(string3(1:(readCount3-2)));
        end
    end    
end

end

