%%%MATLAB plotting code. Makes line plots for various types of data. It
%%%takes in an input variable myVar and saves the graph to outputName. The
%%%string titleString is printed at the top of the graph.

%%%% color blind friendly color palettes when multiple lines are plotted are from 
%%%%https://www.nature.com/articles/nmeth.1618/figures/2

%%%Mostly written by Julia Lazzari-Dean. Thanks
%%%to Ben Adler for font size and axis location printing presets at end of
%%%script.

function [] = makeLinePlot(myVar,outputName,titleString)

% %SIMPLE SINGLE LINE PLOT, no fits
% figure;
% hold on;
% p1 = plot(myVar(:,1),myVar(:,2),'-','Color',[163/255 31/255 1],'Linewidth',0.75);
% 
% %graph formatting
% xlabel('Frame');
% ylabel('Ratio');
% title(titleString);
% %xlim([0 11.5]);
% %Xt = 0:1:13;
% %set(gca,'Xtick',Xt);
% ylim([0.5 2.9]);
% Yt = 0.5:0.5:3;
% set(gca,'Ytick',Yt);

% % %PLOT OF AVERAGED EGF TRACES
% % assuming first column is time, second column is hbss values, third column is
% % hbss sem, fourth column is egf values, fifth column is egf sem
% % for this one, showing both HBSS and EGF on same plot
% figure;
% hold on;
% eb1 = errorbar(myVar(:,1),myVar(:,2),myVar(:,3),'o-');
% eb2 = errorbar(myVar(:,1),myVar(:,4),myVar(:,5),'o-');
% set(eb1,'Color',[0 0 0],'MarkerSize',2,'MarkerFaceColor',[0 0 0],'LineWidth',1);
% set(eb2,'Color',[1 0.2 0],'MarkerSize',2,'MarkerFaceColor',[1 0.2 0],'LineWidth',1);
% 
% % format and labels
% xlabel('Frame Start Time (min)');
% ylabel('Voltage (mV)');
% title('VF2.1.Cl Long Term EGF Responses');
% ylim([-50 -10]);
% xlim([0 15]);
% legend('HBSS','EGF','Location','northeast');
% legend('boxoff');

% EPHYS - PLOT FIRST CELL AND SET UP GRAPH
%plot first cell
myCell=figure;
s1=scatter(myVar(:,1),myVar(:,2),10,'o','filled');
set(s1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
% blue [0 0.445 0.695]
ax = gca;
set(ax.YAxis,'TickLabelFormat','%,.1f');
%adding line of best fit
p = polyfit(myVar(:,1),myVar(:,2),1);
yfit = polyval(p,myVar(:,1));
hold on; 
l1=plot(myVar(:,1),yfit,'LineWidth',1);
set(l1,'Color',[0 0 0]);

%graph formatting
xlabel('Membrane Potential (mV)');
ylabel('Lifetime (ns)');
xlim([-100 100]);
ylim([2.4 3.4]);
set(gca,'YTick',2.4:0.2:3.4);
set(gca,'XTick',-80:40:80);
title(titleString);

% %ADD PATCH NUMBER 2
% s2=scatter(myVar(:,1),myVar(:,3)./1000,10,'o','filled');
% set(s2,'MarkerFaceColor',[0.898 0.621 0],'MarkerEdgeColor',[0.898 0.621 0]);
% p2 = polyfit(myVar(:,1),myVar(:,3),1);
% yfit2 = polyval(p2,myVar(:,1));
% l2=plot(myVar(:,1),yfit2./1000,'LineWidth',1);
% set(l2,'Color',[0.898 0.621 0]);
% %orange [0.898 0.621 0]
% % 
% %ADD PATCH NUMBER 3
% s3=scatter(myVar(:,1),myVar(:,4)./1000,10,'o','filled');
% set(s3,'MarkerFaceColor',[0.844 0.367 0],'MarkerEdgeColor',[0.844 0.367 0]);
% p3 = polyfit(myVar(:,1),myVar(:,4),1);
% yfit3 = polyval(p3,myVar(:,1));
% l3=plot(myVar(:,1),yfit3./1000,'LineWidth',1);
% set(l3,'Color',[0.844 0.367 0]);
% %vermillion [0.844 0.367 0]
% 
% % ADD PATCH NUMBER 4
% s4=scatter(myVar(:,1),myVar(:,5)./1000,10,'o','filled');
% set(s4,'MarkerFaceColor',[0 0.617 0.449],'MarkerEdgeColor',[0 0.617 0.449]);
% p4 = polyfit(myVar(:,1),myVar(:,5),1);
% yfit4 = polyval(p4,myVar(:,1));
% l4=plot(myVar(:,1),yfit4./1000,'LineWidth',1);
% set(l4,'Color',[0 0.617 0.449]);
% %bluish green [0 0.617 0.449]
% 
% % ADD PATCH NUMBER 5
% s5=scatter(myVar(:,1),myVar(:,6)./1000,10,'o','filled');
% set(s5,'MarkerFaceColor',[0.336 0.703 0.910],'MarkerEdgeColor',[0.336 0.703 0.910]);
% p5 = polyfit(myVar(:,1),myVar(:,6),1);
% yfit5 = polyval(p5,myVar(:,1));
% l5=plot(myVar(:,1),yfit5./1000,'LineWidth',1);
% set(l5,'Color',[0.336 0.703 0.910]);
% %sky blue [0.336 0.703 0.910]
% % 
% % ADD PATCH NUMBER 6
% s6=scatter(myVar(:,1),myVar(:,7)./1000,10,'o','filled');
% set(s6,'MarkerFaceColor',[0.797 0.473 0.652],'MarkerEdgeColor',[0.797 0.473 0.652]);
% p6 = polyfit(myVar(1:4,1),myVar(1:4,7),1);
% yfit6 = polyval(p6,myVar(1:4,1));
% l6=plot(myVar(1:4,1),yfit6./1000,'LineWidth',1);
% set(l6,'Color',[0.797 0.473 0.652]);
% %reddish purple [0.797 0.473 0.652]
% % 
% % ADD PATCH NUMBER 7
% s7=scatter(myVar(:,1),myVar(:,8)./1000,10,'o','filled');
% set(s7,'MarkerFaceColor',[0.937 0.891 0.258],'MarkerEdgeColor',[0.937 0.891 0.258]);
% p7 = polyfit(myVar(:,1),myVar(:,8),1);
% yfit7 = polyval(p7,myVar(:,1));
% l7=plot(myVar(:,1),yfit7./1000,'LineWidth',1);
% set(l7,'Color',[0.937 0.891 0.258]);
% %yellow [0.937 0.891 0.258]

% %%I=0 plots
% myCell=figure;
% s1=scatter(myVar(:,1),myVar(:,2),20,'o','filled');
% set(s1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
% xlabel('V_m, Electrophysiology (mV)');
% ylabel('V_m, FLIM (mV)');
% xlim([-45 15]);
% ylim([-45 15]);
% title('A431 Resting Membrane Potential');
% 
% %add a line of best fit
% p = polyfit(myVar(:,1),myVar(:,2),1);
% yfit = polyval(p,myVar(:,1));
% hold on; 
% %just plot from max to min because points are out of order
% plotPoints = zeros(2);
% plotPoints(1,1) = max(myVar(:,1));
% plotPoints(2,1) = min(myVar(:,1));
% indexMax = find(myVar(:,1)==max(myVar(:,1)));
% indexMin = find(myVar(:,1)==min(myVar(:,1)));
% plotPoints(1,2) = yfit(ind2sub(size(myVar(:,1)),indexMax));
% plotPoints(2,2) = yfit(ind2sub(size(myVar(:,1)),indexMin));
% l1=plot(plotPoints(:,1),plotPoints(:,2),':','LineWidth',1);
% set(l1,'Color',[0 0 0]);

%clean up the figure - JACS standard formatting
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

%save to an output file. 
saveas(gcf,[outputName,'.pdf']);

end