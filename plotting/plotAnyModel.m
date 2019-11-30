function [] = plotAnyModel(iStruct,model,oName)

%This function plots the voltage dependence of each of the fit components
%exported into a structure by "reprocessWholeDatasets.m". It currently
%supports exponential decay models with 1-4 terms. Plots are formatted
%according to JACS figure guidelines.

%It also calculates the average line of best fit for the parameter, as well
%as the intra and inter cell RMSD/error values that would accompany using
%that parameter as a proxy for Vm.

if(strcmp(model,'1exp'))
    %list of parameters to process from the dataset
    fitMe = {'tau';'chiSq'};

elseif(strcmp(model,'2exp'))
    %list of parameters to process from the dataset
    fitMe = {'tm';'t1';'a1';'t2';'a2';'chiSq'};    
    
elseif(strcmp(model,'3exp'))
    %list of parameters to process from the dataset
    fitMe = {'tm';'t1';'a1';'t2';'a2';'t3';'a3';'chiSq'};    
    
elseif(strcmp(model,'4exp'))
    %list of parameters to process from the dataset
    fitMe = {'tm';'t1';'a1';'t2';'a2';'t3';'a3';'t4';'a4';'chiSq'};

elseif(strcmp(model,'1ln'))
    %list of parameters to process from the dataset
    fitMe = {'geoMean';'geoStd';'chiSq'};

elseif(strcmp(model,'2ln'))
    error('plotting for this model not yet implemented');
else
    error('Model type not recognized');
end
   
%total number of cells is equal to the number of unique cellIDs
cellIDs = [iStruct(:,1).cellID];
cellIDList = unique(cellIDs)';
nCells = numel(cellIDList);

%partially instantiate the patches structure with the values that will be
%used for all fits
patches = struct('date',"",'patchID',"",'cellID',0,'potentials',0);

%calculate the line of best fit for each parameter vs. Vm, with its corresponding
%slope, yintercept, and r squared.
for i = 1:nCells
    s = iStruct([iStruct.cellID]==cellIDList(i,1)); %identify substructure
    %extract the relevant data into the patches struct
    patches(i,1).potentials = [s(:,1).potential]';
    patches(i,1).chiSq = [s(:,1).chiSq]';
    
    %save the metadata
    patches(i,1).cellID = s(1,1).cellID;
    patches(i,1).patchID = s(1,1).patchID;
    patches(i,1).date = s(1,1).date;
    
    %find the number of saved components
    if (contains(model,'exp'))
        taus = [s(:,1).tau]';
        nExp = numel(taus)/numel(patches(i,1).potentials);
    elseif (contains(model,'ln'))
        nExp = -1;
    else
        error('Unrecognized model string');
    end
    
    %extract the taus and coefficients
    if(nExp == 1)
        patches(i,1).tau = [s(:,1).tau]';
    elseif (nExp == -1)
        if (strcmp(model,'1ln'))
            for j = 1:size(fitMe,1)
                name = fitMe{j,1};
                patches(i,1).(name) = [s(:,1).(name)]';
            end
        elseif (strcmp(model,'1ln'))
            error('LN2 plotting is not yet implemented');
        end
    else
        patches(i,1).tm = [s(:,1).tm]';
        as = [s(:,1).a]';
        
        %extract specific taus
        for j = 1:size(taus,1)/nExp
            for k=1:nExp
                patches(i,1).(['t' num2str(k)])(j,1) = taus(nExp*j-(nExp-k),1);
                patches(i,1).(['a' num2str(k)])(j,1) = as(nExp*j-(nExp-k),1);
            end
        end
    end
    
    %for each of the listed parameters, determine the line of best fit and
    %save its slope, intercept and r squared to the structure
    for j = 1:size(fitMe,1)
        thisVar = fitMe{j,1};
        params = polyfit(patches(i,1).potentials,patches(i,1).(thisVar),1);
        patches(i,1).([thisVar '_m']) = params(1);
        patches(i,1).([thisVar '_b']) = params(2);
        yEval = polyval(params,patches(i,1).potentials);
        r = corrcoef(patches(i,1).(thisVar), yEval);
        patches(i,1).([thisVar '_r2']) = r(1,2)^2;
    end
end

%instantiate up a structure to receive the averages for each parameter
% fitSummary = struct('param',"",'m_avg',0,'m_stdev',0,'m_SEM',0,'b_avg',0,...
%     'b_stdev',0,'b_SEM',0,'r2_avg',0,'r2_stdev',0,'r2_SEM',0,'oRMSD',0,...
%     'ssRMSD_individSlope',0,'ssRMSD_individS_stdev',0,'ssRMSD_individS_SEM',0);
fitSummary = struct('param',"");

%make plots for the listed parameters
for i=1:size(fitMe,1)
    thisVar = fitMe{i,1};
    if (exist('figList','var') == 0)
        figList = figure;
        axList = gca;
    else
        figList(end+1,1) = figure;
        axList(end+1,1) = gca;
    end
    
    %plot the line of best fit for all cells
    for j=1:nCells
        lineVals = [patches(j,1).([thisVar '_m']),patches(j,1).([thisVar '_b'])];
        yEval = polyval(lineVals,patches(j,1).potentials);
        plot(patches(j,1).potentials,yEval,'Color',[0.8 0.8 0.8],'Linewidth',1);
        hold on;
        
        %calculate the intra cell RMSD (using cell-specific slopes, y intercepts)
        ssErrList = 0;
        ssPredV = 0;
        for k = 1:numel(patches(j,1).potentials)
            ssPredV(k,1) = (patches(j,1).(thisVar)(k,1) - patches(j,1).([thisVar '_b']))/patches(j,1).([thisVar '_m']);
            ssErrList(k,1) = ssPredV(k,1) - patches(j,1).potentials(k,1);
        end
        patches(j,1).([thisVar '_ssRMSD']) = sqrt((mean(ssPredV)-mean(patches(j,1).potentials))^2 + var(ssErrList,1));

    end
    
    %aggregate the parameters for the lines of best fit
    fitSummary(i,1).param = thisVar;
    mList = [patches(:,1).([thisVar '_m'])]';
    bList = [patches(:,1).([thisVar '_b'])]';
    r2List = [patches(:,1).([thisVar '_r2'])]';
    ssRMSDList = [patches(:,1).([thisVar '_ssRMSD'])]';    
    %calculate measures of central tendency for these values
    fitSummary(i,1).m_avg = mean(mList);
    fitSummary(i,1).m_stdev = std(mList);
    fitSummary(i,1).m_SEM = std(mList)/sqrt(numel(mList));
    fitSummary(i,1).b_avg = mean(bList);
    fitSummary(i,1).b_stdev = std(bList);
    fitSummary(i,1).b_SEM = std(bList)/sqrt(numel(bList));
    fitSummary(i,1).r2_avg = mean(r2List);
    fitSummary(i,1).r2_stdev = std(r2List);
    fitSummary(i,1).r2_SEM = std(mList)/sqrt(numel(r2List));
    fitSummary(i,1).r2_avg = mean(r2List);
    fitSummary(i,1).r2_stdev = std(r2List);
    fitSummary(i,1).r2_SEM = std(mList)/sqrt(numel(r2List));
    fitSummary(i,1).ssRMSD_avg = mean(ssRMSDList);
    fitSummary(i,1).sSRMSD_stdev = std(ssRMSDList);
    fitSummary(i,1).ssRMSD_SEM = std(ssRMSDList)/sqrt(numel(ssRMSDList));
    
    %calculate the inter cell RMSD
    oErrorList = 0;
    for j = 1:size(mList)
        oErrorList(j,1) = (bList(j,1) - fitSummary(i,1).b_avg)/fitSummary(i,1).m_avg;
    end
    fitSummary(i,1).oRMSD = sqrt(mean(oErrorList)^2 + var(oErrorList,1));

    %plot the a line with the average slope and intercept
    lineVals = [fitSummary(i,1).m_avg, fitSummary(i,1).b_avg];
    yEval = polyval(lineVals,-80:40:80);
    plot(-80:40:80,yEval,'Color',[0 0 0],'Linewidth',1);
    
    %add titles, adjust ranges
    title([oName '_' thisVar],'Interpreter','none');
    xlabel('Potential (mV)','HorizontalAlignment','center');
    xlim([-100 100]);
    
    %adjust y axis based on what is being plotted
    bAvg = round(fitSummary(i,1).b_avg,1);
    if (contains(thisVar,'geo') || contains(thisVar,'t'))
        if(fitSummary(i,1).b_avg < 0.5)
            yll = 0;
            yul = 1;
        else
            yll = bAvg - 0.5;
            yul = bAvg + 0.5;
        end
        ylim([yll yul]);
        yt = yll:0.2:yul;
        set(gca,'YTick',yt);
        ylabel([thisVar ' (ns)'],'HorizontalAlignment','center');
    elseif (contains(thisVar,'a'))
        if(fitSummary(i,1).b_avg < 0.2)
            yll = 0;
            yul = 0.4;
        else
            yll = bAvg - 0.2;
            yul = bAvg + 0.2;
        end
        ylim([yll yul]);
        yt = yll:0.1:yul;
        set(gca,'YTick',yt);
        ylabel(thisVar,'HorizontalAlignment','center');
    else
        ylabel(thisVar,'HorizontalAlignment','center');        
    end    
    %will let chi squared auto scale because that will vary a lot
end

for i=1:size(figList,1)
    %generally clean up the plot
    set(get(axList(i,1),'xlabel'),'FontSize', 10,'Color',[0 0 0]);
    set(get(axList(i,1),'ylabel'),'FontSize', 10,'Color',[0 0 0]);
    set(get(axList(i,1),'title'),'FontSize', 10,'Color',[0 0 0],'FontWeight','normal');
    set(axList(i,1),'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
    set(axList(i,1),'TickDir','in')
    set(axList(i,1),'box','off');
    set(axList(i,1),'LineWidth',1);
    set(axList(i,1),'FontSize',8);
    set(figList(i,1),'color','w');
    
    %Format Positioning On the Page
    pbaspect(axList(i,1), [1 1 1]);
    axList(i,1).Units = 'inches';
    axList(i,1).Position = [0.4 0.4 1.75 1.75];
    
    %save each figure as pdf
    name = axList(i,1).Title.String;
    saveas(figList(i,1),[name,'.pdf']);
end

%save the patches structure
save(oName,'patches');
save(oName,'fitSummary','-append');

end

