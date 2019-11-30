%%comment this function here

function [] = plotESMCombo(data,outputName)
%set color scheme
zeroColor = [0.5 0.5 0.5];
negEightyColor = [205/256 0 205/256];

%set up a figure with subplots
figure;
xlabel('Tau (ns)');
ylabel('Coefficient Weight');
hold on;
[patchIDs, ~] = unique([data.patchID].', 'rows', 'stable');
nPatches = length(patchIDs);
subplot(nPatches,1,1);
axList = gca;

for i=1:size(patchIDs,1)
    %separate out the data from this patchID
    count = 0;
    for j = 1:size(data,1)
        if(data(j,1).patchID == patchIDs(i,1))
            if (count == 0)
                thisData = data(j,1);
            else
                thisData(end+1,1) = data(j,1);
            end
            count = count+1;
        end
    end
       
    %put in dashed or straight lines based on whether this was recorded
    %first or second
    
    %plot the ESM result
    subplot(nPatches,1,i);
    axList(end+1,1) = gca;
    hold on;
    title(strcat(thisData(1,1).date,"  ",thisData(1,1).patch));
    for j=1:size(thisData,1)
        %determine the line color based on the potential
        if(thisData(j,1).potential == 0)
            color = zeroColor;
        elseif (thisData(j,1).potential == -80)
            color = negEightyColor;
        else
            err('Unrecognized potential');
        end
        plot(0.05:0.05:6,thisData(j,1).aFinal,'-','Color',color,'LineWidth',0.75);
        ylabel('Weight');
        xlim([0 6]);
        ylim([0 1]);
    end
    
    %add a legend if this was the first subplot
    if i==1 && thisData(1,1).potential == 0
        legend({"0 mV","-80 mV"});
        legend('boxoff');
    elseif i==1 && thisData(1,1).potential == -80
        legend({"-80 mV","0 mV"});    
        legend('boxoff');
    end
end

%clean up the axes
for i=1:size(axList,1)
    set(get(axList(i,1),'xlabel'),'FontSize', 10,'Color',[0 0 0]);
    set(get(axList(i,1),'ylabel'),'FontSize', 10,'Color',[0 0 0]);
    set(get(axList(i,1),'title'),'FontSize', 10,'Color',[0 0 0]);
    set(axList(i,1),'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
    set(axList(i,1),'TickDir','in')
    set(axList(i,1),'LineWidth',0.75);
    set(axList(i,1),'FontSize',10);
    set(axList(i,1),'YTick',[0 0.25 0.5 0.75 1]);
    pbaspect(axList(i,1),[6 1 1]);
end

%add an overall x label to the bottom graph
xlabel('Tau (ns)');
sgtitle(outputName,'Interpreter','None','Color',[0 0 0]);

%Format Positioning On the Page
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [7 7]);
set(gcf,'PaperPosition',[1 0.5 6.5 5.5]);
set(gcf,'PaperPositionMode','Manual');
set(gcf,'color','w');
box off;

%save to an output file. switch which line is commented to switch styles.
saveas(gcf,[outputName,'.pdf']);


end

