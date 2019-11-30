%This function takes a full dataset from various days and processes it
%using batchFloptimizeGeneral.

%The input structure should have the following fields:
    %1. Date
    %2. IRF
    %3. array of decays to be fit with that IRF.
    %4. potentials that correspond to each decay
    %5. list of the number of potentials for each patch (related to the
    %number of patches)

%In this iteration, it is agnostic to the potential/patch for each, just
%reducing the total number of files/variables/function calls.

function [] = reprocessWholeDatasets(iStruct,model,mode,oName,cShift,shiftFixed,dcOffset)

if(strcmp(mode,'ephys'))
    %separate out the input structure such that it can be added to the output
    %structures and metadata retained. in doing this, also check the input
    %structure and validate that the number of decays, patches, and potentials
    %checks out.
    md = struct('date',"",'cellID',"","patchID",0,'potential',0);
    count = 1;
    cellNum = 1;
    
    %determine the total number of cells and assign a unique cell number to
    %each patch
    for i=1:size(iStruct,1)
        IDarray = cellNum:1:(cellNum - 1 + numel(iStruct(i,1).potPerP));
        iStruct(i,1).cellIDs = IDarray';
        cellNum = cellNum + numel(iStruct(i,1).potPerP);
    end
    nCells = cellNum - 1;
    
    for i = 1:size(iStruct,1)
        nDecays = sum(iStruct(i,1).potPerP);
        iStruct(i,1).nDecays = nDecays;
        nPotentials = numel(iStruct(i,1).potentials);
        nCells = numel(iStruct(i,1).potPerP);
        nIDs = numel(iStruct(i,1).patchIDs);
        if (nCells ~= nIDs)
            error(['Entry ' num2str(i) ' has the wrong number of patch IDs.']);
        elseif (nPotentials ~= nDecays)
            error(['Entry ' num2str(i) ' has the wrong number of potentials listed.']);
        end
        
        for j=1:nDecays
            md(count,1).date = iStruct(i,1).date;
            md(count,1).potential = iStruct(i,1).potentials(j,1);
            for k = 1:nIDs
                if(j <= sum(iStruct(i,1).potPerP(1:k,1)))
                    md(count,1).patchID = iStruct(i,1).patchIDs(k,1);
                    md(count,1).cellID = iStruct(i,1).cellIDs(k,1);
                    break;
                end
            end
            count = count + 1;
        end
    end
elseif(strcmp('single',mode))
    %do i need to do anything here? i think it's actually fine since each
    %row has one decay
    md = struct('date',"");
    for i=1:size(iStruct,1)
        nDecays = size(iStruct(i,1).decays,2);
        if nDecays > 1
            disp('Warning - more than one decay per image in single mode.');
        end
        iStruct(i,1).nDecays = nDecays;
    end
end

%if model = all, fit with all of the options and save out all the data
if (strcmp(model,'all'))
    %process everything with a 1 exponential fit
    for i = 1:size(iStruct,1)
        if i == 1
            exp1fit = batchFloptimizeGeneral(iStruct(i,1).decays,iStruct(i,1).IRF,...
                '1exp',cShift,shiftFixed,dcOffset,2);
        else
            nDecays = iStruct(i,1).nDecays;
            result = batchFloptimizeGeneral(iStruct(i,1).decays,...
                iStruct(i,1).IRF,'1exp',cShift,shiftFixed,dcOffset,2);
            exp1fit(end+1:end+nDecays,1) = result;
        end
    end
    
    %process everything with a 2 exponential fit
    for i = 1:size(iStruct,1)
        if i == 1
            exp2fit = batchFloptimizeGeneral(iStruct(i,1).decays,iStruct(i,1).IRF,...
                '2exp',cShift,shiftFixed,dcOffset,[1 1 0.9 2.5],[0 0 0 0]);
        else
            nDecays = iStruct(i,1).nDecays;
            result = batchFloptimizeGeneral(iStruct(i,1).decays,...
                iStruct(i,1).IRF,'2exp',cShift,shiftFixed,dcOffset,...
                [1 1 0.9 2.5],[0 0 0 0]);
            exp2fit(end+1:end+nDecays,1) = result;
        end
    end
    
    %process everything with a 3 exponential fit
    for i = 1:size(iStruct,1)
        if i == 1
            exp3fit = batchFloptimizeGeneral(iStruct(i,1).decays,iStruct(i,1).IRF,...
                '3exp',cShift,shiftFixed,dcOffset,[1 1 1 0.1 0.9 2.5],[0 0 0 0 0 0]);
        else
            nDecays = iStruct(i,1).nDecays;
            result = batchFloptimizeGeneral(iStruct(i,1).decays,...
                iStruct(i,1).IRF,'3exp',cShift,shiftFixed,dcOffset,...
                [1 1 1 0.1 0.9 2.5],[0 0 0 0 0 0]);
            exp3fit(end+1:end+nDecays,1) = result;
        end
    end
    
    %process everything with a 4 exponential fit
    for i = 1:size(iStruct,1)
        if i == 1
            exp4fit = batchFloptimizeGeneral(iStruct(i,1).decays,iStruct(i,1).IRF,...
                '4exp',cShift,shiftFixed,dcOffset,[1 1 1 1 0.1 0.5 1 3],[0 0 0 0 0 0 0 0]);
        else
            nDecays = iStruct(i,1).nDecays;
            result = batchFloptimizeGeneral(iStruct(i,1).decays,...
                iStruct(i,1).IRF,'4exp',cShift,shiftFixed,dcOffset,...
                [1 1 1 1 0.1 0.5 1 3],[0 0 0 0 0 0 0 0]);
            exp4fit(end+1:end+nDecays,1) = result;
        end
    end
    
    %process everything with a LN dist fit
    for i = 1:size(iStruct,1)
        if i == 1
            ln1fit = batchFloptimizeGeneral(iStruct(i,1).decays,iStruct(i,1).IRF,...
                'LogNorm1',cShift,shiftFixed,dcOffset,[1 0.1],[0 0]);
        else
            nDecays = iStruct(i,1).nDecays;
            result = batchFloptimizeGeneral(iStruct(i,1).decays,...
                iStruct(i,1).IRF,'LogNorm1',cShift,shiftFixed,dcOffset,...
                [1 0.1],[0 0]);
            ln1fit(end+1:end+nDecays,1) = result;
        end
    end
    
    %%CURRENTLY THROWING ERRORS
%     %process everything with a sum of 2 LN dist fit
%     for i = 1:size(iStruct,1)
%         if i == 1
%             ln2fit = batchFloptimizeGeneral(iStruct(i,1).decays,iStruct(i,1).IRF,...
%                 'LogNorm2',cShift,shiftFixed,dcOffset,[1 1 0.5 1 0.1 0.1],[0 0 0 0 0 0]);
%         else
%             nDecays = sum(iStruct(i,1).potPerP);
%             result = batchFloptimizeGeneral(iStruct(i,1).decays,...
%                 iStruct(i,1).IRF,'LogNorm2',cShift,shiftFixed,dcOffset,...
%                 [1 1 0.5 1 0.1 0.1],[0 0 0 0 0 0]);
%             ln2fit(end+1:end+nDecays,1) = result;
%         end
%     end
    
end

%add the metadata back to each output structure - workaround by converting
%to tables and merging tables
% Convert structures to tables
mdT = struct2table(md);
exp1fitT = struct2table(exp1fit);
exp2fitT = struct2table(exp2fit);
exp3fitT = struct2table(exp3fit);
exp4fitT = struct2table(exp4fit);
ln1fitT = struct2table(ln1fit);

% Concatenate tables
merge1exp = [mdT,exp1fitT];
merge2exp = [mdT,exp2fitT];
merge3exp = [mdT,exp3fitT];
merge4exp = [mdT,exp4fitT];
merge1ln = [mdT,ln1fitT];

% Convert table to structure
exp1results = table2struct(merge1exp);
exp2results = table2struct(merge2exp);
exp3results = table2struct(merge3exp);
exp4results = table2struct(merge4exp);
ln1results = table2struct(merge1ln);

%write all of these tables to csv
writetable(removevars(merge1exp,{'decay','irf','residTrace'}),[oName '_1exp.csv']);
writetable(removevars(merge2exp,{'decay','irf','residTrace'}),[oName '_2exp.csv']);
writetable(removevars(merge3exp,{'decay','irf','residTrace'}),[oName '_3exp.csv']);
writetable(removevars(merge4exp,{'decay','irf','residTrace'}),[oName '_4exp.csv']);
writetable(removevars(merge1ln,{'decay','irf','residTrace'}),[oName '_1ln.csv']);

%write all of these tables to csv in addition to saving the variables
%save variables to output file
save(oName,'exp1results');
save(oName,'exp2results', '-append');
save(oName,'exp3results', '-append');
save(oName,'exp4results', '-append');
save(oName,'ln1results', '-append');

% %also do this for the double log normal data
% ln2fitT = struct2table(ln2fit);
% mergeln2 = [mdT,ln2fitT];
% ln2results = table2struct(mergeln2);
% writetable(ln2results,[oName '_2ln'],'writevariablenames',0);
% save(oName,'ln2results', '-append');

end

