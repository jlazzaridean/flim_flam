%This function fits an exponential decay to each pixel of a photon histogram imported with importData.m

%The .bin or .sdt photon histograms must have already been parsed and
%connected to the relevant IRFs and metadata by importData. The function
%will render certain fit results as an image and save the corresponding
%TIFF file in the same path specified within the input structure. This
%function supports binning, both in the spatial dimension and in nanosecond
%time. The function returns a structure with the fit results at each pixel
%and also saves this structure to a .mat file. 

%Assessing fit quality: Two measures of fit quality
%are generated graphically - [1] histograms of the reduced chi squared at each
%pixel for the first image from each cID and [2] a plot of the weighted
%residuals for the first image from each cID (averaged across all pixels).
%These fit quality graphs are saved to the same folder as the original
%data.

%Written by Julia Lazzari-Dean. Last edits October 12, 2019.

%sample function calls:
    %myPxFits = pxwiseAnalysis(myImportedData,["tau";"chiSq";"overlay_tau"],...
    %  'imageSaveFolderName',2,'std',1,chang1expc0,[2.5 3.0]);
    
    %myPxFits = pxwiseAnalysis(myImportedData,["tau";"photons_raw";"photons"],...
    %   'imageSaveFolderName2',1,'bh',2,chang1expc0);
    
%Descripton of input parameters:
    %data - the output of 'pxwise' mode run inputData
    %render - parameters for which to generate TIFF images, as an N x 1
    %list of strings. The images will be saved to the newly generated 
    %directory oFolder within the path name saved with the input data structure.
        %Images are generally saved as an X by Y array of the value. If an
        %overlay with the photon image is required, for example with tm (weighted average lifetime),
        %you would use the string "overlay_tm" Overlays can be generated
        %for any parameter. For each overlay, you will need to specify a
        %range as an additional parameter in varargin (see below).
        %To render the binned photons image, use "photons". To render the
        %unbinned photons image, use "photons_raw". To render the sum of
        %squared errors, use "SSE." To render the reduced chi squared, use
        %"chiSq."
    %oFolder - the directory created within the path name in the data
    %structure where images will be saved
    %spatialBin - the extent of binning. The meaning of this depends on the
    %binType (following parameter)
    %binType - either 'bh' or 'std', specifies how the binning is completed
        %'bh' - binning similar to that in SPCImage. bin 1 = data from all
        %pixels one away are added to a given pixel (photon upsampling;
        %lifetime is effectively a moving average)
        %'std' - a more conventional binning scheme. bin 2 would indicate
        %that a 2x2 chunk of pixels is added together to make 1 binned
        %pixel
    %adcBin - the binning in the ADC dimension of the data. The binning
    %format on this is always 'std'
    %config - an imported structure for a config file. This can either be
    %generated directly by the helper function readConfig() or used from
    %the output of the importData function.
    %varargin - additional parameters are required whenever an overlay is
    %requested. For the example above, to display a tm overlay on the
    %photon image with tm scaled from 1 to 2 ns, varargin would be [1 2].
    %If you requested overlay_tm from 1-2 ns and overlay_tau1 from 0.3 to
    %0.8 ns, varargin would have two range sets and be written [1 2], [0.3 0.8]
    %in whichever order you specified the strings in render.


%Note: All overlays will be rendered over the photon count image, with the LUTs
%optimized for the range of each image. To render over a specific photon
%count range, call the helper function renderImage directly and specify the
%parameter photonMax. For 'std,' the lifetime is overlaid on the binned
%photon image. In 'bh' binning, the upsampled lifetime map is overlaid on
%the unbinned photon image (as in SPCImage).

%Dependencies (all *.m): writeTIFF, renderImages, renderOverlay,
%binRawTCSPC, convol, assignParams_2exp, assignParams_3exp,
%assignParams_4exp, floptimize3_1exp, floptimize3_2exp, floptimze3_3exp,
%floptimize4_3exp, sortATs

function [data] = pxwiseAnalysis(data,render,oFolder,spatialBin,binType,adcBin,config,varargin)

%add subfolders to the path
addpath('parsingTracing');
addpath('fittingVariousModels');

%varargin will be the range for each overlay as [limLow limHigh] in the
%order they are listed in the render matrix
if(~isempty(render))
    nOverlays = nnz(contains(render,'overlay'));
    if(nOverlays ~= size(varargin,2))
        error('Please enter a range as [low high] for each overlay requested in varargin');
    end
end

%parse the data in the config file to get the necessary parameters
model = config(1,1).model;
period_ns = config(1,1).period_ns;
cShift = config(1,1).cShift;
shiftFixed = config(1,1).shiftFixed;
startParam = config(1,1).startParam;
fixedParam = config(1,1).fixedParam;
thresh = config(1,1).threshold;
stFi = config(1,1).stFi./adcBin;
offFixed = config(1,1).offFixed;
ADCres = config(1,1).ADCres/adcBin;
viewDecay = config(1,1).viewDecay;

%identify the number of exponential components to use in generating a
%structure of the correct size
nExp = str2double(model(1));
nExp = round(nExp);
fsep = filesep;

%create the directory oFolder in the path mentioned by the data
folderInfo = dir(data(1,1).pName);
if(folderInfo(1,1).isdir)
    outPath = strcat(data(1,1).pName,fsep,oFolder);
else
    disp('Invalid directory. Saving to current matlab path.')
    outPath = oFolder;
end
mkdir(outPath);

%iterate over the number of images to fit
for i=1:size(data,1)
    data(i,1).model = model; %save the model for the record
    if ((spatialBin > 1 || adcBin > 1) && strcmp(binType,'std'))
        [binnedImg, fitIRF] = binRawTCSPC(data(i,1).tcspc,data(i,1).IRF,spatialBin,binType,adcBin);
    elseif ((spatialBin > 0 || adcBin > 1) && strcmp(binType,'bh'))
        [binnedImg, fitIRF] = binRawTCSPC(data(i,1).tcspc,data(i,1).IRF,spatialBin,binType,adcBin);
    else
        binnedImg = data(i,1).tcspc;
        fitIRF = data(i,1).IRF;
    end
    %create a photons matrix and write this to data
    data(i,1).photonsBinned = sum(binnedImg,3);
    data(i,1).binType = binType;
    data(i,1).spatialBin = spatialBin;
    data(i,1).adcBin = adcBin;
    fitMask = data(i,1).photonsBinned > thresh; %create a mask of indices to fit
    data(i,1).photons = sum(data(i,1).tcspc,3);    
    imgP = permute(binnedImg,[3 1 2]); %permute the binnedImg to make slicing faster
    
    %go across the image and fit all pixels that are above threshold
    for j=1:size(imgP,2)
        for k=1:size(imgP,3)
            %call the appropriate function
            if fitMask(j,k)  == 1
                %fitMe = reshape(binnedImg(j,k,:), [numel(binnedImg(j,k,:)) 1 1]);
                fitMe = imgP(:,j,k);
                switch model
                    case '1exp'
                        [tF, cF, offset, chiSq, residTrace, SSE, ...
                            exitFlag] = floptimize3_1exp(fitMe,fitIRF,...
                            cShift,shiftFixed,offFixed,startParam,stFi,period_ns,viewDecay);
                    case '2exp'
                        [tm, aFs, tFs, cF, offset, chiSq, residTrace, SSE, ...
                            exitFlag] = floptimize3_2exp(fitMe,fitIRF,...
                            cShift,shiftFixed,offFixed,startParam,fixedParam,stFi,period_ns,viewDecay);
                    case '3exp'
                        [tm, aFs, tFs, cF, offset, chiSq, residTrace, SSE, ...
                            exitFlag] = floptimize3_3exp(fitMe,fitIRF,...
                            cShift,shiftFixed,offFixed,startParam,fixedParam,stFi,period_ns,viewDecay);
                    case '4exp'
                        [tm, aFs, tFs, cF, offset, chiSq, residTrace, SSE, ...
                            exitFlag] = floptimize3_4exp(fitMe,fitIRF,...
                            cShift,shiftFixed,offFixed,startParam,fixedParam,stFi,period_ns,viewDecay);
                    otherwise
                        error('Unrecognized model string.');
                end
                if(exitFlag < 1) %the solver did not converge
                    error(['Exit Flag ' num2str(exitFlag) ' for fit. Exiting.']);
                end
                %add these values to the correct arrays
                if(strcmp(model,'1exp'))
                    data(i,1).tau(j,k) = tF;
                elseif(strcmp(model,'2exp') || strcmp(model,'3exp') || strcmp(model,'4exp'))
                    data(i,1).tau(j,k,1:nExp) = reshape(tFs,[1 1 nExp]);
                    data(i,1).a(j,k,1:nExp) = reshape(aFs,[1 1 nExp]);
                    data(i,1).tm(j,k) = tm;
                end
                data(i,1).SSE(j,k) = SSE;
                data(i,1).chiSq(j,k) = chiSq;
                data(i,1).cShift(j,k) = cF;
                data(i,1).offset(j,k) = offset;
                data(i,1).resid(j,k,:) = reshape(residTrace,[1 1 numel(residTrace)]);
            else
                if(strcmp(model,'1exp'))
                    data(i,1).tau(j,k) = NaN;
                elseif(strcmp(model,'2exp') || strcmp(model,'3exp') || strcmp(model,'4exp'))
                    data(i,1).tau(j,k,1:nExp) = NaN;
                    data(i,1).a(j,k,1:nExp) = NaN;
                    data(i,1).tm(j,k) = NaN;
                end
                data(i,1).SSE(j,k) = NaN;
                data(i,1).chiSq(j,k) = NaN;
                data(i,1).cShift(j,k) = NaN;
                data(i,1).offset(j,k) = NaN;
                data(i,1).resid(j,k,1:ADCres) = NaN;
            end
        end
    end
    
    %render the images - one at a time in this case so that the user can
    %see as things unfold
    if nOverlays > 0
        renderImages(data(i,1),render,outPath,-1,varargin{1,:});
    else
         renderImages(data(i,1),render,outPath,-1);
    end
    
    close all
end

%save the result as a .mat file
oName = strcat(outPath,fsep,data(1,1).date,'_analyzedData');
save(oName,'data','-v7.3');

%close all open figures
close all

%show the weighted chi squared pixelwise histogram for the first image for each cID (this is
%imperfect but will give a rough sense of what is going on) - break this up
%into sets of 6 plots per exported image
cIDList = [data(:,1).cID]';
cIDUnique = unique(cIDList);
cIDCount = numel(cIDUnique);
nFigs = ceil(cIDCount/6);
%instantiate a gobject with 2*nFigs
figList = gobjects(nFigs*2,1);
axList = gobjects(cIDCount*2,1);
plotCount = 1;
%iterate through the required number of figures and do the plotting
for i=1:nFigs
    figList(i,1) = figure;
    %get ax for each subplot
    for j=1:6
        if(plotCount > cIDCount) %on the last loop, stop before you get to 6
            break;
        end
        axList(plotCount,1) =  subplot(3,2,j);
        %take the subset of the data with this CID
        subset = data([data.cID]' == cIDUnique(plotCount,1));
        %plot the first decay in this subset
        thisChiSq = reshape(subset(i,1).chiSq,[numel(subset(i,1).chiSq) 1]);
        thisChiSq (isnan(thisChiSq)) = [];
        histogram(thisChiSq);
        xlabel('Reduced Chi Squared');
        ylabel('Abundance');
        title(['cID ' num2str(cIDUnique(plotCount,1)) ', 1st Image']);
        plotCount = plotCount + 1;
    end
end

%also show the average residuals (is this even a reasonable thing?) for the
%first image of each coverslip
plotCount = 1;
%calculate the time resolution
period_ns = config(1,1).period_ns;
time = 0:period_ns/ADCres:(period_ns - period_ns/ADCres);
time = time';
%iterate through the required number of figures and do the plotting
for i=1:nFigs
    figList(i+nFigs,1) = figure;
    %get ax for each subplot
    for j=1:6
        if(plotCount > cIDCount) %on the last loop, stop before you get to 6
            break;
        end
        axList(cIDCount + plotCount,1) = subplot(3,2,j);
        %take the subset of the data with this CID
        subset = data([data.cID]' == cIDUnique(plotCount,1));
        %plot the first decay in this subset
        avgResid = nanmean(nanmean(subset(1,1).resid,1),2);
        avgResid = reshape(avgResid,[numel(avgResid) 1 1]);
        plot(time,avgResid,'LineWidth',1);
        xlabel('Time (ns)');
        ylabel('Avg. Wtd. Res.');
        title(['cID ' num2str(cIDUnique(plotCount,1)) ', 1st Image']);
        plotCount = plotCount + 1;
    end
end

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
    pbaspect(axList(i,1), [2.5 1 1]);
end

names = ["_chiSqHist";"_avgResid"];
fs = filesep;
fullOName = strcat(outPath,fs,data(i,1).date,'_pxWise');
for i=1:size(figList,1)
    set(figList(i,1),'color','w');
    if (i <= size(figList,1)/2)
        saveas(figList(i,1),strcat(fullOName,names(1,1),num2str(i),'.pdf'));
    else
        saveas(figList(i,1),strcat(fullOName,names(2,1),num2str(i-nFigs),'.pdf'));        
    end
end

end

