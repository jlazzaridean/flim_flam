function [bData] = frameBinDecays(data,fBin)
%this function takes in a structure with data, retains some metadata (just
%copies the first entry - so only one time series per structure for now)
%and adds the number of decays together specified by fBin
% fieldsToCopy = {"cID";"imageID";"cellType";"conc_nM";"irf_index";
%     "fName";"pName";"IRF";"irfName";"date"};
fieldsToCopy = {"condition";"cID";"imageID";"cellType";"conc_nM";"irf_index";
    "fName";"pName";"IRF";"irfName";"date"};

%identify the unique cIDs and imageIDs in the structure (only want to bin a
%recording; don't want to bin the whole structure if it contains multiple
%recordings)
recordings = struct("cID",0,"imageID",0,"nFrames",0);
cIDList = [data.cID]';
cIDUnique = unique(cIDList);
count = 1;
for j = 1:size(cIDUnique,1)
    cIDMatch = data([data.cID]' == cIDUnique(j,1));
    iIDList = [cIDMatch.imageID]';
    iIDUnique = unique(iIDList);
    for k = 1:size(iIDUnique,1)
        %slice and get only these frames
        iIDMatch = cIDMatch([cIDMatch.imageID]' == iIDUnique(k,1));
        %populate the structure with this information
        recordings(count,1).cID = cIDUnique(j,1);
        recordings(count,1).imageID = iIDUnique(k,1);
        recordings(count,1).nFrames = size(iIDMatch,1);
        count = count + 1;
    end
end

framesList = [recordings.nFrames]';
bFramesList = floor(framesList/fBin);
nbFrames = sum(bFramesList,1);
bData(1:nbFrames,1) = struct('frameID',0,'decays',0);
countRec = 1;
countFrames = 1;

for j = 1:size(cIDUnique,1)
    cIDMatch = data([data.cID]' == cIDUnique(j,1));
    iIDList = [cIDMatch.imageID]';
    iIDUnique = unique(iIDList);
    for k = 1:size(iIDUnique)
        iIDMatch = cIDMatch([cIDMatch.imageID]' == iIDUnique(k,1));
        bFrames = bFramesList(countRec,1);
        %iterate over binned frames and generate summed decays
        for i=1:bFrames
            summedDecay = zeros(size(iIDMatch(1,1).decays));
            for m = 0:(fBin-1)
                summedDecay = summedDecay + iIDMatch(i*fBin+m-fBin+1,1).decays;
            end
            bData(countFrames,1).decays = summedDecay;
            bData(countFrames,1).frameID = i;
            for m = 1:size(fieldsToCopy)
                bData(countFrames,1).(fieldsToCopy{m,1}) = iIDMatch(1,1).(fieldsToCopy{m,1});
            end
            countFrames = countFrames + 1;
        end
        
        countRec = countRec + 1;
    end
end

end

