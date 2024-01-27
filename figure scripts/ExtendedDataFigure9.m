
%% A,B
% Plot the pplace preference Vectors of tracked feature coders for all transitions
DoPlot = [1 0];
load(['B1-B2.mat'])
numses = [1 2];
helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);
[plcVentral1, plcDorsal1,DNrestActFract1,pooledf(1)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,0);

load(['B3-T1.mat'])
numses = [1 2];
helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);
[plcVentral2, plcDorsal2,DNrestActFract2,pooledf(2)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,0);

load([pathname 'B3-Tn.mat'])
numses = [1 2];
helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);
[plcVentral3, plcDorsal3,DNrestActFract3,pooledf(3)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,0);

load([pathname 'Tn-P.mat'])
numses = [1 2];
helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);
[plcVentral4, plcDorsal4,DNrestActFract4,pooledf(4)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,0);

%% C
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[15 5 [18.3 10]],...
    'PaperUnits','centimeters',...
    'PaperSize', [18.3 10])

subplot(1,2,1)
binSize = 5;
bins = -75-binSize/2:binSize:75+binSize/2;

plcVcount1 = histcounts(plcVentral1,bins,'normalization','probability');
plcVcount2 = histcounts(plcVentral2,bins,'normalization','probability');
plcVcount3 = histcounts(plcVentral3,bins,'normalization','probability');
plcVcount4 = histcounts(plcVentral4,bins,'normalization','probability');

hold off
x = bins(2:end);%+binSize/2;
plot([0 0],[0 1],'--','color',[.5 .5 .5])
hold on
a = plot(x,cumsum(plcVcount1));
a(1) = plot(x,cumsum(plcVcount1));
a(2) = plot(x,cumsum(plcVcount2));
a(3) = plot(x,cumsum(plcVcount3));
a(4) = plot(x,cumsum(plcVcount4));
% plot(x,normcdf(x,0,10),'r-')


xlim([x(1) x(end)])
ylim([0 1])
xlabel('PF shift (cm)')
ylabel('cum. prob.')
legend(a, {'B1-B2' 'B3-T1' 'T1-Tn' 'Tn-P1'},'location','southeast')
legend('boxoff')
title('PF shift Ventral')
% grid on
box off

subplot(1,2,2)

plcDcount1 = histcounts(plcDorsal1,bins,'normalization','probability');
plcDcount2 = histcounts(plcDorsal2,bins,'normalization','probability');
plcDcount3 = histcounts(plcDorsal3,bins,'normalization','probability');
plcDcount4 = histcounts(plcDorsal4,bins,'normalization','probability');

hold off
x = bins(2:end);%+binSize/2;
plot([0 0],[0 1],'--','color',[.5 .5 .5])
hold on
a = plot(x,cumsum(plcDcount1));
a(1) = plot(x,cumsum(plcDcount1));
a(2) = plot(x,cumsum(plcDcount2));
a(3) = plot(x,cumsum(plcDcount3));
a(4) = plot(x,cumsum(plcDcount4));
% plot(x,normcdf(x,0,19),'r-')

xlim([x(1) x(end)])
ylim([0 1])
xlabel('PF shift (cm)')
ylabel('cum. prob.')
legend(a, {'B1-B2' 'B3-T1' 'T1-Tn' 'Tn-P1'},'location','southeast')
legend('boxoff')
title('PF shift Dorsal')
% grid on
box off
