%This function fits exponential decays to a structure of global decays imported by importData.

%For each entry in the structure, the code will fit exponential decays as
%described in the input configuration structure (all fit model parameters
%must be specified there as documented in the manual). If more than one
%decay is present for each entry (e.g. if BTM ROIs were traced on the
%image), each decay will be fit. A structure with the fit results is
%returned.

%The script will save the results to a .mat file specified by the parameter
%oName. The script will also save an 'unrolled' version of the data to
%.csv, where metadata is duplicated such that each row in the .csv
%corresponds to the results of fitting an individual decay. 

%This function generates two graphical measures of fit quality. First, a
%histogram of all the chi squared values for all of the fits performed.
%Second, a plot of the weighted residuals for all of the decays for the
%first entry for each cID. Both of these figures are saved into the same
%folder as the other results.

%Written by Julia Lazzari-Dean. Last edits October 12, 2019.

%Sample function call:
    % globallyFitData = fitDataStruct(importedData,'myOutputName',chang1expc0);

%Input parameters:
    %data - the output of importData.m from any of the global modes
    %('single','global_btm' or 'global_fiji')
    %oName - the base output name for saving results. These files will be
    %saved to the path specified in the data structure (where the .sdt and
    %.bin files are located)
    %configS - a structure containing the configuration and system
    %parameters for the data. It may be used as an output from either
    %readConfig or importData.
    
%Dependencies (all .m): batchFloptimizeExp, sortATs, assignParams_2exp,
%assignParams_3exp, assignParams_4exp, convol, floptimize3_1exp,
%floptimize3_2exp, floptimize3_3exp, floptimize3_4exp

function [data] = fitDataStruct(data,oName,configS)

%add subfolders to the path
addpath('fittingVariousModels');

%parse the data in the config file to get the necessary parameters
model = configS(1,1).model;
modelChar = char(model);
nExp = str2double(modelChar(1,1));

%iterate through the data and fit
for i = 1:size(data,1)
    %save relevant metadata about the fit process to the output structure
    data(i,1).model = model;
    data(i,1).fixedParam = configS(1,1).fixedParam;
    data(i,1).startParam = configS(1,1).startParam;
    if(isempty(data(i,1).decays))
        error('No decays found for index %d',i);
    end
    data(i,1).totalPhotons = sum(data(i,1).decays,1)';
    
    %bin the IRF and the decay if the user wanted
    data(i,1).adcBin = configS(1,1).adcBin;
    [binnedDecays,binnedIRF] = binADCGlobal(data(i,1).decays, data(i,1).IRF,data(i,1).adcBin);
    data(i,1).bDecays = binnedDecays;
    data(i,1).bIRF = binnedIRF;
    %call the fit functions and add the results back to the 
    theseResults = batchFloptimizeExp(binnedDecays,binnedIRF,configS);
    theseFields = fieldnames(theseResults);
    %write the as and the taus in their own loop
    for j = 1:size(theseFields,1)
        thisField = theseFields{j,1};
        if(strcmp(thisField,'tau'))
            for k = 1:size(theseResults,1)
                data(i,1).tau(k,1:nExp) = theseResults(k,1).tau(1,1:nExp);
            end
        elseif(strcmp(thisField,'a'))
            for k = 1:size(theseResults,1)
                data(i,1).a(k,1:nExp) = theseResults(k,1).a(1,1:nExp);
            end
        else
            data(i,1).(thisField) = [theseResults(:,1).(thisField)]';
        end
    end
end

%expand the structure for writing to .csv
fieldsAll = fieldnames(data);
if(any(strcmp(fieldsAll,'masks'))) %data came from images
    dataToWrite = rmfield(data,{'bIRF','bDecays','residTrace','irfName','masks'});
else
    dataToWrite = rmfield(data,{'bIRF','bDecays','residTrace','irfName'});
end

%remove taus, as, start and fixed param from multiexp models - will write them later
count = 1;
if(~strcmp(configS(1,1).model,'1exp')) 
    fieldToWrite = fieldnames(rmfield(dataToWrite,{'tau','a','startParam','fixedParam'}));
else
    fieldToWrite = fieldnames(dataToWrite);
end

expandedData = struct('date',"");
for i=1:size(dataToWrite,1)
    nFits = size(dataToWrite(i,1).chiSq,1); %each row is a fit. different columns of tau are different components.
    for j=1:nFits
        for k=1:size(fieldToWrite,1)
            thisF = fieldToWrite{k,1};
            if(size(dataToWrite(i,1).(thisF),1) == 1)
                expandedData(count,1).(thisF) = dataToWrite(i,1).(thisF);
            else
                expandedData(count,1).(thisF) = dataToWrite(i,1).(thisF)(j,1);
            end
        end
        %if it is a >1exp model, write the taus and a's separately
        if(~strcmp(configS(1,1).model,'1exp'))
            expandedData(count,1).tau(1,1:nExp) = dataToWrite(i,1).tau(j,1:nExp);
            expandedData(count,1).a(1,1:nExp) = dataToWrite(i,1).a(j,1:nExp);
            expandedData(count,1).startParam = configS(1,1).startParam;
            expandedData(count,1).fixedParam = configS(1,1).fixedParam;
        end
        count = count + 1;
    end
end

%% FIT QUALITY METRICS

%first, figure out how many figures we will need
%then, show the weighed residuals for the first decay for each cID (this is
%imperfect but will give a rough sense of what is going on) - break this up
%into sets of 6 plots per exported image
cIDList = [expandedData(:,1).cID]';
cIDUnique = unique(cIDList);
cIDCount = numel(cIDUnique);
nFigs = ceil(cIDCount/6);
figList = gobjects(nFigs+1,1);
axList = gobjects(cIDCount+1,1);
%show the user some fit quality metrics in graphical format
%first, take a histogram of all of the chi squared values in the
%expandedData structure
chiSqAll = [expandedData(:,1).chiSq]';
figList(1,1) = figure;
axList(1,1) = gca;
histogram(chiSqAll);
title('All Chi Squared Values for Data');
xlabel('Reduced Chi Squared');
ylabel('Decay Count');

%now make the weighted residuals plots
plotCount = 1;
%calculate the time resolution
adcRes = configS(1,1).adcRes/configS(1,1).adcBin;
nsPeriod = configS(1,1).nsPeriod;
time = 0:nsPeriod/adcRes:(nsPeriod - nsPeriod/adcRes);
%iterate through the required number of figures and do the plotting
for i=1:nFigs
    figList(i+1,1) = figure;
    %get ax for each subplot
    for j=1:6
        if(plotCount > cIDCount) %on the last loop, stop before you get to 6
            break;
        end
        axList(plotCount+1,1) = subplot(3,2,j);
        %take the subset of the data with this CID
        subset = data([data.cID]' == cIDUnique(plotCount,1));
        %plot the first decay in this subset
        plot(time,subset(1,1).residTrace,'LineWidth',1);
        xlabel('Time(ns)');
        xlim([0 nsPeriod]);
        ylabel('Weighted Residual');
        title(['cID ' num2str(cIDUnique(plotCount,1)) ', 1st Entry']);
        plotCount = plotCount + 1;
    end
end

%% SAVING FIGURES AND DATA
%save data in the path of the first entry of the structure
thisPath = data(1,1).pName;
fsep = filesep;
if(thisPath(end) ~= fsep)
    thisPath(end+1) = fsep;
end
fullOName = strcat(thisPath,oName);

%clean up and save the figures
for i=1:size(axList,1)
    %generally clean up the plot
    set(get(axList(i,1),'xlabel'),'FontSize', 10,'Color',[0 0 0]);
    set(get(axList(i,1),'ylabel'),'FontSize', 10,'Color',[0 0 0]);
    set(get(axList(i,1),'title'),'FontSize', 10,'Color',[0 0 0],'FontWeight','normal');
    set(axList(i,1),'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
    set(axList(i,1),'TickDir','in')
    set(axList(i,1),'box','off');
    set(axList(i,1),'LineWidth',1);
    set(axList(i,1),'FontSize',8);
    if i > 1
        pbaspect(axList(i,1), [2.5 1 1]);
    end
end

%Format Positioning On the Page - only for the first plot
pbaspect(axList(1,1), [1 1 1]);
axList(1,1).Units = 'inches';
axList(1,1).Position = [0.4 0.4 1.75 1.75];

names = ["_chiSqHist";"_sampleResid"];
for i=1:size(figList,1)
    set(figList(i,1),'color','w');
    if i==1
        saveas(figList(i,1),strcat(fullOName,names(1,1),'.pdf'));
    else
        saveas(figList(i,1),strcat(fullOName,names(2,1),num2str(i-1),'.pdf'));        
    end
end

%write the reshaped data to .csv (for easy pandas export - one line per fit
%result)
writetable(struct2table(expandedData),[fullOName '_fit.csv']);

%also save the unexpanded structure to .mat
save([fullOName '_fit'],'data','configS','-v7.3');

end

