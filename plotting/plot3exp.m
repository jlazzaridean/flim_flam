function [] = plot3exp(iStruct,oName)

%this function plots the voltage dependence of each individual parameter of
%a 3 component fit (as exported into a structure by the
%"reprocessWholeDatasets.m" function).

%list of parameters to process from the datasets (editing this here only
%may cause issues later down)
fitMe = {'tm';'t1';'a1';'t2';'a2';'t3';'a3';'chiSq'};

%total number of cells is equal to the number of unique cellIDs
cellIDs = [iStruct(:,1).cellID];
cellIDList = unique(cellIDs)';
nCells = numel(cellIDList);

%set up a structure to put the results for each cell
patches = struct('date',"",'patchID',"",'cellID',0,'potentials',0,'tm',0,'t1',0,...
    'a1',0,'t2',0,'a2',0,'t3',0,'a3',0,'chiSq',0,'tm_m',0,'tm_b',0,'tm_r2',0,...
    't1_m',0,'t1_b',0,'t1_r2',0,'a1_m',0,'a1_b',0,'a1_r2',0,'t2_m',0,'t2_b',0,...
    't2_r2',0,'a2_m',0,'a2_b',0,'a2_r2',0,'t3_m',0,'t3_b',0,'t3_r2',0,'a3_m',0,...
    'a3_b',0,'a3_r2',0,'chiSq_m',0,'chiSq_b',0,'chiSq_r2',0);

%calculate the line of best fit for each parameter vs. Vm, with its corresponding
%slope, yintercept, and r squared.
for i = 1:nCells
    s = iStruct([iStruct.cellID]==cellIDList(i,1)); %identify substructure
    %extract the relevant data into the patches struct
    patches(i,1).potentials = [s(:,1).potential]';
    patches(i,1).tm = [s(:,1).tm]';
    patches(i,1).chiSq = [s(:,1).chiSq]';
    
    %save the metadata
    patches(i,1).cellID = s(1,1).cellID;
    patches(i,1).patchID = s(1,1).patchID;
    patches(i,1).date = s(1,1).date;
    
    %extract the taus and coefficients
    taus = [s(:,1).tau]';
    as = [s(:,1).a]';
    for j = 1:size(taus,1)/3
        if (patches(i,1).t1 == 0)
            patches(i,1).t1(1,1) = taus(3*j-2,1);
            patches(i,1).t2(1,1) = taus(3*j-1,1);
            patches(i,1).t3(1,1) = taus(3*j,1);
            patches(i,1).a1(1,1) = as(3*j-2,1);
            patches(i,1).a2(1,1) = as(3*j-1,1);
            patches(i,1).a3(1,1) = as(3*j,1);
        else
            patches(i,1).t1(end+1,1) = taus(3*j-2,1);
            patches(i,1).t2(end+1,1) = taus(3*j-1,1);
            patches(i,1).t3(end+1,1) = taus(3*j,1);
            patches(i,1).a1(end+1,1) = as(3*j-2,1);
            patches(i,1).a2(end+1,1) = as(3*j-1,1);
            patches(i,1).a3(end+1,1) = as(3*j,1);
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
    end
    
    %plot the average line of best fit
    
    %add titles, adjust ranges
    title([oName '_' thisVar],'Interpreter','none');
    xlabel('Potential (mV)','HorizontalAlignment','center');
    ylabel(thisVar,'HorizontalAlignment','center');
    xlim([-100 100]);

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
    set(figList(i,1),'PaperUnits','inches');
    set(figList(i,1),'PaperSize', [2.5 2.5]);
    set(figList(i,1),'PaperPosition',[0 0 2.5 2.5]);
    set(figList(i,1),'PaperPositionMode','Manual');
    
    %save each figure as pdf
    name = axList(i,1).Title.String;
    saveas(figList(i,1),[name,'.pdf']);
end

%save/determine averages for each of the parameters?

%save the patches structure
save(oName,'patches');


end

