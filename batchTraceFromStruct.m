%This function allows averaging lifetimes across regions of interest in data that has already been fit pixel by pixel
%by the function pxwiseAnalysis.m. Results are saved to .mat format and also to .csv in the
%file path indicated within the input structure.

%Written by Julia Lazzari-Dean. Last edits October 12, 2019. Content
%borrows heavily from my previously developed batchTraceMembranes2 suite.

%Sample function call:
    %batchTraceFromStruct(myPxWiseFitData,'tau','conc_nM',["cellType";"incubationTime"],...
    %  'myOutputFolderName','myOutputFileName');

%Input parameter descriptions:
    %data - This is the structure output generated from 'pxwise' mode of
    %importData followed by processing with pxwiseAnalysis.m. The data must
    %contain the fields cellType and smallestObject for the batch tracing
    %to work properly, as well as the field binType to render the correct
    %version of the photons image.
    
    %signalParam - the field of the fit data that you would like to use as
    %an output. A common selection would be "tau" for a monoexponential fit
    %or "tm" for a multiexponential fit. If you would like to access a
    %particular coefficient or component of a multiexponential fit, refer
    %to it as, for example, "a1" or "a2" or "tau2" (don't refer to
    %components outside the range of what exists in your fit model).
    
    %categoryParam - the metadata field for grouping results. This can be
    %any field that you supplied in the original metadata file.
    
    %importantParam - metadata fields that you would like to retain in the
    %output .csv for processing in other software. For instance, if all of
    %your data are described by two categories (e.g. treatment
    %and incubationTime), these values will also be saved with the average
    %values of signalParam and the corresponding categoryParam for each
    %generated ROI. These fields should be input as an N by 1 array of
    %strings.
    
    %outputFolderName - a folder name to create within the path named in the
    %data structure.
    
    %outputFileName - the base filename for summaries of the results (box
    %plots of the data by category and by coverslip, .csv outputs of the
    %data for each ROI).
    
%Overall, this function finds cells in photon count
%images using Otsu's method of thresholding. It then identifies isolated
%groups of cells using Matlab's built-in connectivity assessment
%(bwconncomp). From these connectivity maps, it generates regions of
%interest (ROIs) for analysis of resting membrane potential in the
%corresponding lifetime image. The user is asked to approve ROIs generated for
%every image (rainbow color scheme delineates cell groups) by
%pressing "Enter" when the image appears. If any debris hasn't been
%filtered out (i.e. if debris is color-coded), it can be removed by a left
%click on the offending ROI before the image is approved with "Enter."    

%The script then applies these photon count image ROIs to signalParam image. 
%The mean value of lifetime in each ROI is tabulated for each image,
%coverslip, and category (entered in "metadata", see below). These values
%are saved to .mat files in an output directory, and they are also exported
%to .csv for input into Pandas dataframes (python) or other convenient
%plotting software.

%The script generates two box plots, one of the data broken down by
%coverslip, and the other of the data broken down by category. These plots
%are also saved to the output directory within the path specified by the 
%input data strucure. Binary TIFF images of all of the masks generated are
%also saved to a masks subfolder within the input data structure.

%Dependencies (all *.m): ccToMask, findROIsToMerge, maskImage2_otsu,
%traceMembranes4, writeTIFF

function [] = batchTraceFromStruct(data,signalParam,categoryParam,importantParam,outputFolderName,outputFileName)

%add subfolders to the path
addpath('parsingTracing');

%the data should contain smallest_object and cell_type (from the initial
%metadata)
fields = fieldnames(data);
if(~any(strcmp(fields,'smallestObject')))
    error('Input pxwise data structure should contain smallest_object metadata');
elseif(~any(strcmp(fields,'cellType')))
    error('Input pxwise data structure should contain cell_type metadata');
end

%make an output directory and a directory within that for the masks
fsep = filesep;
outPath=strcat(data(1,1).pName,fsep,outputFolderName);
mkdir(outPath);
maskPath = strcat(outPath,fsep,'masks');
mkdir(maskPath);

%total n samples
measCount = 1;

%open the files, trace them, and save the results
for i = 1:size(data,1)
    %call the tracing script until user approves the ROIs
    if(strcmp(data(i,1).binType,'std'))
        goodROIs = 0;
        while(~goodROIs)
            [masks, goodROIs] = traceMembranes4(data(i,1).photonsBinned,data(i,1).cellType,data(i,1).smallestObject);
        end
    elseif(strcmp(data(i,1).binType,'bh'))
        goodROIs = 0;
        while(~goodROIs)
            [masks, goodROIs] = traceMembranes4(data(i,1).photons,data(i,1).cellType,data(i,1).smallestObject);
        end
    else
        error('binType field in the data must be either bh or std');
    end
    
    %save the masks to the structure and also as BW images
    %apply the masks to the data and save the results as meanSignal for
    %each ROI in the main data structure
    data(i,1).signalParam = signalParam;
    data(i,1).masks = masks;
    data(i,1).meanSignal = 0;
    oFName = strcat(data(i,1).fName(1:end-4),'_mask');
    if(strcmp(signalParam,'tau1') || strcmp(signalParam,'tau2') || strcmp(signalParam,'tau3') || strcmp(signalParam,'tau4'))
        sigChar = char(signalParam);
        tNum = str2double(sigChar(end));
        signalImg = data(i,1).tau(:,:,tNum);
    elseif(strcmp(signalParam,'a1') || strcmp(signalParam,'a2') || strcmp(signalParam,'a3') || strcmp(signalParam,'a4'))
        sigChar = char(signalParam);
        aNum = str2double(sigChar(end));
        signalImg = data(i,1).a(:,:,aNum);        
    else
        signalImg = data(i,1).(signalParam);
    end
    for j=1:size(masks,3)
        writeTIFF(data(i,1).masks(:,:,j),'analysis',maskPath,strcat(oFName,num2str(j)));
        image = immultiply(masks(:,:,j),signalImg);
        image(image == 0) = NaN;
        data(i,1).meanSignal(j,1) = nanmean(nanmean(image,1),2);
        measCount = measCount + 1;
    end
           
    close all;
end

%convert category to strings before transpose to avoid shenanigans
if(ischar(data(1,1).(categoryParam)))
    for i=1:size(data,1)
        charArray = data(i,1).(categoryParam);
        data(i,1).(categoryParam) = convertCharsToStrings(charArray);
    end
end

%find the unique categories and coverslip IDs listed in the metadata cell array
categories = [data.(categoryParam)]';
uniqueCategories = unique(categories);
cIDs = [data.cID]';
uniqueCIDs = unique(cIDs);

%instantiate another structure that will hold the summarized results with
%less metadata, as well as some summary structures
if(isa(uniqueCategories(1,1),'string'))
    results(1:measCount-1,1) = struct('cID',0,'fName',"",signalParam,0,categoryParam,"");
    categoryStruct(1:size(uniqueCategories,1),1) = struct(categoryParam,"",'mean',0,'stdev',0,'sem',0,'count',0);
else
    results(1:measCount-1,1) = struct('cID',0,'fName',"",signalParam,0,categoryParam,0);
    categoryStruct(1:size(uniqueCategories,1),1) = struct(categoryParam,0,'mean',0,'stdev',0,'sem',0,'count',0);
end
plotCategory = zeros(1,size(uniqueCategories,1));
plotCoverslip = zeros(1,size(uniqueCIDs,1));
coverslipStruct(1:size(uniqueCIDs,1),1) = struct('cID',0,'mean',0,'stdev',0,'sem',0,'count',0);

%unpack the results into another structure where one value is in each line
writeCount = 1;
for i=1:size(data,1)
    for j=1:size(data(i,1).meanSignal,1)
        results(writeCount,1).(signalParam) = data(i,1).meanSignal(j,1);
        results(writeCount,1).cID = data(i,1).cID;
        for k=1:size(importantParam,1)
            results(writeCount,1).(importantParam(k,1)) = data(i,1).(importantParam(k,1));
        end
        results(writeCount,1).fName = data(i,1).fName;
        results(writeCount,1).(categoryParam) = data(i,1).(categoryParam);
        writeCount = writeCount + 1;
    end
end

%generate summary statistics and an array to plot for each category
for i = size(uniqueCategories,1):-1:1
    if(isa(uniqueCategories(i,1),'string'))
        subset = results(strcmp([results.(categoryParam)]',uniqueCategories(i,1)));
    else
        subset = results([results.(categoryParam)]' == uniqueCategories(i,1));
    end
    categoryStruct(i,1).(categoryParam) = uniqueCategories(i,1);
    categoryStruct(i,1).mean = nanmean([subset(:,1).(signalParam)]);
    categoryStruct(i,1).stdev = nanstd([subset(:,1).(signalParam)]);
    categoryStruct(i,1).count = numel([subset(:,1).(signalParam)]);
    categoryStruct(i,1).sem = categoryStruct(i,1).stdev/sqrt(categoryStruct(i,1).count);
    plotCategory(1:categoryStruct(i,1).count,i) = [subset(:,1).(signalParam)]';
end

%generate summary statistics and an array to plot for each cID
for i = size(uniqueCIDs,1):-1:1
    subset = results([results.cID]' == uniqueCIDs(i,1));
    coverslipStruct(i,1).cID = uniqueCIDs(i,1);
    coverslipStruct(i,1).mean = nanmean([subset(:,1).(signalParam)]);
    coverslipStruct(i,1).stdev = nanstd([subset(:,1).(signalParam)]);
    coverslipStruct(i,1).count = numel([subset(:,1).(signalParam)]);
    coverslipStruct(i,1).sem = coverslipStruct(i,1).stdev/sqrt(coverslipStruct(i,1).count);
    plotCoverslip(1:coverslipStruct(i,1).count,i) = [subset(:,1).(signalParam)]';
end

%display these values in a table
disp(strcat("Summary statistics for ",signalParam," in category ",categoryParam));
disp(struct2table(coverslipStruct));
disp(struct2table(categoryStruct));

%remove zeros from the plot arrays
plotCoverslip(plotCoverslip == 0) = NaN;
plotCategory(plotCategory == 0) = NaN;

%graph the results by coverslip
f1 = figure;
figList = f1;
ax1 = gca;
axList = ax1;

%make the boxplot by coverslip
boxplot(plotCoverslip);
xlabel('Coverslip ID');
ylabel(signalParam,'Interpreter','none');
title([outputFileName ' - by Coverslip'],'interpreter','none'); %this name-value pair disables underscores as subscripts
xticklabels(uniqueCIDs);

%print the sample sizes for each coverslip on the graph
for i = 1:size(coverslipStruct,1)
    szTxt = ['n = ' num2str(coverslipStruct(i,1).count)];
    text(i,1100,szTxt,'HorizontalAlignment','center');
end

%graph the results by category
f2 = figure;
figList(end+1) = f2;
ax2 = gca;
axList(end+1) = ax2;

%make the boxplot by category
boxplot(plotCategory);
xlabel(categoryParam,'Interpreter','none');
ylabel(signalParam,'Interpreter','none');
title([outputFileName ' - by Category'],'interpreter','none'); %this name-value pair disables underscores as subscripts
xticklabels(uniqueCategories);

%print the sample sizes for each condition on the graph
for i = 1:size(categoryStruct,1)
    szTxt = ['n = ' num2str(categoryStruct(i,1).count)];
    text(i,1100,szTxt,'HorizontalAlignment','center');
end


%clean up the figures
for i=1:size(figList,2)
    set(figList(i),'PaperUnits','inches');
    set(figList(i),'PaperSize', [7 6]);
    set(figList(i),'PaperPosition',[1 0.5 6.5 5.5]);
    set(figList(i),'PaperPositionMode','Manual');
    set(figList(i),'color','w');
end

for i=1:size(axList,2)
    set(get(axList(i),'xlabel'),'FontSize', 14,'Color',[0 0 0]);
    set(get(axList(i),'ylabel'),'FontSize', 14,'Color',[0 0 0]);
    set(get(axList(i),'title'),'FontSize', 16,'Color',[0 0 0]);
    set(axList(i),'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
    set(axList(i),'TickDir','in')
    set(axList(i),'box','off');
    set(axList(i),'LineWidth',0.75);
    set(axList(i),'FontSize',14);
end

%save the plots
overallOutputName = strcat(outPath,fsep,outputFileName);
f1Name = strcat(overallOutputName,'_byCoverslip');
saveas(f1,[char(f1Name),'.pdf']);
f2Name = strcat(overallOutputName,'_byCat');
saveas(f2,[char(f2Name),'.pdf']);

%save the results as .mat files and the unrolled results structure to .csv
save(overallOutputName, 'data','-v7.3');
save(overallOutputName, 'coverslipStruct', '-append');
save(overallOutputName, 'categoryStruct', '-append');
writetable(struct2table(results),strcat(overallOutputName,'_results.csv'));

end