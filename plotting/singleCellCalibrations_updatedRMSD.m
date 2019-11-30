%This script takes a bunch of single cell voltage step calibration data,
%fits them with linear fits, and records those fit parameters. It then
%generates a plot with the linear fit for each cell in gray. It also calculates the average by potential (as
%well as some other stats) and plots this in black on top of the graph.

% ***This script is currently hardcoded to work with potentials -80, -40, 0
%,40, and 80. It will include other potentials in the fit for the
%individual cells, but it will exclude them from any sort of 'by potential'
%averaging to make the overall line of best fit.

%RMSD calculations used in VF-FLIM: overall (inter cell) is oRMSD, 'tracking' (intra cell)
%is ssRMSD_indS_avg (SEM reported for this is ssRMSD_indS_SEM).

%One of the main advantages of this script is that it does not require any
%calculations to be done in Excel in advance. All of the calculations from
%the voltage steps spreadsheets are made and saved in this script.

%This script takes in three input variables. dataList is a list of the
%potential and lifetime pairs from a bunch of different cells. The first
%column is potentials and the second column is lifetimes. All recordings
%for a cell are grouped together, and they are separated from the next cell
%by a zero in both columns (potential and lifetime). Note that the script
%distinguishes between two different cells by looking for a zero in the
%lifetime column. All input data should be in ps, and all output data will
%be in ns. All voltage throughout the script is in mV.

%The other two input variables are more trivial - gTitle is the title for
%the graphs the function will generate, and output is the name of the file
%that all of the analysis/variables/statistics will be saved to. The
%function returns a list of the voltage error (in mV) at each of the
%voltage individual voltage steps submitted. (Note: I tried to do an
%average per cell as well, but that ends up being almost exactly zero
%always by the definition of linear regression.)

%Determinations of Error: The script also estimates the between cell (oRMSD) and
%within cell (ssRMSD) voltage errors. For the between cell error, the
%script determines the RMSD for the y-intercept of each cell as '0 mV'
%relative to the overall calibration for that cell type.
%For the within cell errors, the script calculates the RMSD between the
%applied step and the calibration for that specific cell. For
%ssError_ssAvg, the script gets an RMSD per cell and then averages all of
%those. ssError_ssSEM is the SEM of these values.
%For ssError_overall, the script calculates the RMSD using each
%voltage step as an independent measurement (so all pooled at once, instead
%of treating each cell as an independent measurement).
%All RMSD calculations include both the bias and the variance.

%written by Julia Lazzari-Dean 
%last edited May 23, 2019 - fixing RMSD calculations. Also doing RMSD
%calculations with the average slope (new kind of intra cell error) - now
%differentiate between _indS (individual slope) and _avgS (average
%slope)intra-cell errors.

%features to add: 
%1.trends in single cell error at a potential (to evaluate how well LJP is equilibrating)

function [] = singleCellCalibrations_updatedRMSD(dataList,gTitle,output) 

%parse the input data
%it will be coming in as the same format as an Excel file
i = 1;
cCount = 0;
patches = struct('data',0,'slope',0,'yinter',0,'rsquared',0,'errors_indS',0,...
    'errors_avgS',0,'predV_indS_intra',0,'predV_avgS_intra',0,'oError',0,...
    'ssRMSD_indS',0,'ssRMSD_avgS',0);
overallList = -1;
while (i <= size(dataList,1))
    cCount = cCount + 1;
    mCount = 1;
    while (i <= size(dataList,1) && dataList(i,2) ~= 0)
        patches(cCount).data(mCount,1) = dataList(i,1);
        patches(cCount).data(mCount,2) = dataList(i,2);
        if(overallList(1,1) == -1)
            overallList(1,1) = dataList(i,1);
            overallList(1,2) = dataList(i,2);
        else
            overallList(end+1,1) = dataList(i,1);
            overallList(end,2) = dataList(i,2);
        end
        mCount = mCount + 1;
        i = i + 1;
    end
    i = i + 1;
end

%set up an initial figure
f1 = figure;
ax1 = gca;
hold on;
title(gTitle);
ylabel('Lifetime (ns)');
xlabel('Membrane Potential (mV)');
xlim([-100 100]);
ylim([2.4 3.4]);
set(gca,'XTick',-80:40:80);
set(gca,'YTick',2.4:0.2:3.4);
set(ax1.YAxis,'TickLabelFormat','%,.1f');

%for each entry, calculate a line of best fit and save those parameters
%also plot all of these individual lines in gray onto the best fit plot
for i = 1:size(patches,2)
    p = polyfit(patches(i).data(:,1),patches(i).data(:,2),1);
    patches(i).slope = p(1);
    patches(i).yinter = p(2);
    yfit = polyval(p,patches(i).data(:,1));
    r = corrcoef(patches(i).data(:,2), yfit);
    patches(i).rsquared = r(1,2)^2;
    plot(patches(i).data(:,1),yfit,'Color',[0.8 0.8 0.8],'Linewidth',0.75); 
end

%apply the single-cell lines of best fit and determine errors. save these
%both to the struct and to the overall error list
for i = 1:size(patches,2)
    for j = 1:size(patches(i).data,1)
        patches(i).predV_indS_intra(j) = (patches(i).data(j,2)-patches(i).yinter)/patches(i).slope;
        patches(i).errors_indS(j) = patches(i).predV_indS_intra(j) - patches(i).data(j,1);
    end
    %this RMSD determination includes bias (which should be zero) and
    %variance.
    patches(i).ssRMSD_indS = sqrt((mean(patches(i).predV_indS_intra(:))-mean(patches(i).data(:,1)))^2 + var(patches(i).errors_indS(:),1));
end
%determine the RMSD using each measurement independently
ssRMSD_indS_overall = sqrt((mean([patches.predV_indS_intra])-mean(overallList(:,1)))^2 + var([patches.errors_indS],1));
ssRMSD_indS_avg = mean([patches.ssRMSD_indS]);
ssRMSD_indS_SEM = std([patches.ssRMSD_indS])/sqrt(numel([patches.ssRMSD_indS]));

%sort the overall data by potential so stats can be determined
overallFit = struct('potential',-1,'data',[],'average',-1,'median',-1,'stdev',-1,'sem',-1,'count',-1);
for i = 1:5
   overallFit(i).potential = i*40 - 120; 
end
for i=1:size(overallList,1)
    voltage = overallList(i,1);
    if (mod(voltage,40) ~=0)
        disp('Potential that is not a multiple of 40 was ignored.');
        continue;
    end
    if (isempty(overallFit([overallFit.potential] == voltage).data))
        overallFit([overallFit.potential] == voltage).data(1,1) = overallList(i,2);
    else
        overallFit([overallFit.potential] == voltage).data(end+1,1) = overallList(i,2);
    end
end

%calculate statistics at each potential
for i = 1:5
    overallFit(i).average = round(nanmean(overallFit(i).data),3);
    overallFit(i).median = round(nanmedian(overallFit(i).data),3);
    overallFit(i).stdev = round(nanstd(overallFit(i).data),5);
    overallFit(i).count = numel(overallFit(i).data)-sum(isnan(overallFit(i).data));
    overallFit(i).sem = round(overallFit(i).stdev/sqrt(overallFit(i).count),5);
end

%calculate the line of best fit through the averages at each potential, and
%plot this on top of the existing figure.
%(this script is basically a better version of ephysManyLines)
numPotentials = 5;
averagesForFit = zeros(numPotentials,1);
errorsForFit = zeros(numPotentials,1);
for i = 1:numPotentials
    averagesForFit(i,1) = overallFit(i).average;
    errorsForFit(i,1) = overallFit(i).sem;
end
p = polyfit([-80;-40;0;40;80],averagesForFit(:,1),1);
yfit = polyval(p,[-80;-40;0;40;80]);
eb1 = errorbar([-80;-40;0;40;80], yfit, errorsForFit(:,1),'-o');
set(eb1,'LineWidth',1,'MarkerSize',2,'Color',[0 0 0],'MarkerFaceColor',[0 0 0],'CapSize',4);
r = corrcoef(averagesForFit(:,1), yfit);
rsquaredOverall = r(1,2)^2;
slopeOverall = p(1);
yinterOverall = p(2);

%use the overall slope to determine an RMSD resulting only from the between
%cell y-intercept variability (this is the best estimate of expected error
%in calculating voltage changes in time series)
avgM = mean([patches.slope]);
avgB = mean([patches.yinter]);
semM = std([patches.slope])/sqrt(mCount);
semB = std([patches.yinter])/sqrt(mCount);
for i = 1:size(patches,2)
    for j = 1:size(patches(i).data,1)
        patches(i).predV_avgS_intra(j) = (patches(i).data(j,2)-patches(i).yinter)/avgM;
        patches(i).errors_avgS(j) = patches(i).predV_avgS_intra(j) - patches(i).data(j,1);
    end
    %this RMSD determination includes bias (which should be zero) and
    %variance.
    patches(i).ssRMSD_avgS = sqrt((mean(patches(i).predV_avgS_intra(:))-mean(patches(i).data(:,1)))^2 + var(patches(i).errors_avgS(:),1));
end
%determine the RMSD using each measurement independently
ssRMSD_avgS_overall = sqrt((mean([patches.errors_avgS])-mean(overallList(:,1)))^2 + var([patches.errors_avgS],1));
ssRMSD_avgS_avg = mean([patches.ssRMSD_avgS]);
ssRMSD_avgS_SEM = std([patches.ssRMSD_avgS])/sqrt(numel([patches.ssRMSD_avgS]));

%use the overall calibration to calculate the overall errors from the 0 mV
%points.
for i = 1:size(patches,2)
    %voltageFromInter should be '0' for a perfect calibration; i.e. the
    %yintercept of the single cell should also give a zero with the
    %yintercept of the overall calibration
    patches(i).oError = (patches(i).yinter - yinterOverall)/slopeOverall;
end
oRMSD = sqrt(mean([patches.oError])^2 + var([patches.oError],1));

%clean up the figure
set(get(gca,'xlabel'),'FontSize', 10,'Color',[0 0 0]);
set(get(gca,'ylabel'),'FontSize', 10,'Color',[0 0 0]);
set(get(gca,'title'),'FontSize', 10,'Color',[0 0 0],'FontWeight','normal');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
set(gca,'TickDir','out','TickLength',[0.02 0.035])
set(gca,'box','off');
set(gca,'LineWidth',1);
set(gca,'FontSize',8);
set(gca,'Units','inches','Position',[0.4 0.4 1.4 1.4])
set(gcf,'color','w');

%save data to output file
save(output, 'patches','avgM','avgB','semM','semB');
save(output,'dataList', '-append');
save(output,'overallList', '-append');
save(output,'cCount', '-append');
save(output,'slopeOverall', '-append');
save(output,'rsquaredOverall', '-append');
save(output,'yinterOverall', '-append');
save(output,'overallFit', '-append');
save(output,'ssRMSD_indS_overall','-append');
save(output,'ssRMSD_indS_avg','-append');
save(output,'ssRMSD_indS_SEM','-append');
save(output,'ssRMSD_avgS_overall','-append');
save(output,'ssRMSD_avgS_avg','-append');
save(output,'ssRMSD_avgS_SEM','-append');
save(output,'oRMSD','-append');

outputName1 = strcat(output,'_linesPlot');
saveas(f1,[outputName1,'.pdf']);

end