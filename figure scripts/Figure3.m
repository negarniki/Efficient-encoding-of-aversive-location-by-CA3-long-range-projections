% Figure 3

experiment = {'B1-B2', 'B3-T1', 'B3-Tn', 'Tn-P1'}';
%% Pool preference vectors of tracked components
DoPlot = [0 0];
load(['B1-B2.mat'])
numses = [1 2];
PV = helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);

load(['B3-T1.mat'])
numses = [1 2];
PV(:,:,2) = helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);

load([pathname 'B3-Tn.mat'])
numses = [1 2];
PV(:,:,3) = helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);

load([pathname 'Tn-P.mat'])
numses = [1 2];
PV(:,:,4) = helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);

%% B,F
DoPlot = [0 1];
load([pathname 'Tn-P.mat'])
numses = [1 2];
helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);

%% D,H
DoPlot = [1 0];
load([pathname 'B3-Tn.mat'])
numses = [1 2];
helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);
load([pathname 'Tn-P.mat'])
numses = [1 2];
helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);
%% C
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[3 5 [15 15]],...
    'PaperUnits','centimeters',...
    'PaperSize', [15 15])

x = 1:4;
plot(x,permute(PV(7,2,:),[3 1 2]),'r-o')
hold on
plot(x,permute(PV(7,4,:),[3 1 2]),'k-*')
xlim([0 5])
ylim([0 1])
legend({'I-D' 'D-D'},'location','northwest')
legend("boxoff")
ylabel('Fraction of axons')
ax = gca;
ax.XTick = x;
ax.XTickLabel = experiment;
box off

%% E
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[3 5 [15 15]],...
    'PaperUnits','centimeters',...
    'PaperSize', [15 15])

x = 1:4;
errorbar(x,permute(PV(1,2,:),[3 1 2]),permute(PV(3,2,:),[3 1 2]),'r-o')
hold on
errorbar(x,permute(PV(1,4,:),[3 1 2]),permute(PV(3,4,:),[3 1 2]),'k-*')
xlim([0 5])
ylim([0 1])
legend({'I-D' 'D-D'},'location','northwest')
legend("boxoff")
ylabel('Average Vector length')
ax = gca;
ax.XTick = x;
ax.XTickLabel = experiment;
box off

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[3 5 [15 15]],...
    'PaperUnits','centimeters',...
    'PaperSize', [15 15])

x = 1:4;
plot(x,permute(PV(4,2,:),[3 1 2]),'r-o')
hold on
plot(x,permute(PV(4,4,:),[3 1 2]),'k-*')
xlim([0 5])
ylim([-180 180])
legend({'I-D' 'D-D'},'location','northeast')
legend("boxoff")
ylabel('Average angular distance')
ax = gca;
ax.XTick = x;
ax.XTickLabel = experiment;
box off

%% G
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[3 5 [15 15]],...
    'PaperUnits','centimeters',...
    'PaperSize', [15 15])

x = 1:4;
plot(1:4,permute(PV(7,1,:),[3 1 2]),'r-o')
hold on
plot(1:4,permute(PV(7,3,:),[3 1 2]),'k-*')
xlim([0 5])
ylim([0 1])
legend({'I-D' 'D-D'},'location','northwest')
legend("boxoff")
ylabel('Fraction of axons')
ax = gca;
ax.XTick = x;
ax.XTickLabel = experiment;
box off

%% I
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[3 5 [15 15]],...
    'PaperUnits','centimeters',...
    'PaperSize', [15 15])

x = 1:4;
errorbar(x,permute(PV(1,1,:),[3 1 2]),permute(PV(3,1,:),[3 1 2]),'r-o')
hold on
errorbar(x,permute(PV(1,3,:),[3 1 2]),permute(PV(3,3,:),[3 1 2]),'k-*')
xlim([0 5])
ylim([0 1])
legend({'I-D' 'D-D'},'location','northwest')
legend("boxoff")
ylabel('Average Vector length')
ax = gca;
ax.XTick = x;
ax.XTickLabel = experiment;
box off

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[3 5 [15 15]],...
    'PaperUnits','centimeters',...
    'PaperSize', [15 15])

x = 1:4;
plot(x,permute(PV(4,1,:),[3 1 2]),'r-o')
hold on
plot(x,permute(PV(4,3,:),[3 1 2]),'k-*')
xlim([0 5])
ylim([-180 180])
legend({'I-D' 'D-D'},'location','northeast')
legend("boxoff")
ylabel('Average angular distance')
ax = gca;
ax.XTick = x;
ax.XTickLabel = experiment;
box off