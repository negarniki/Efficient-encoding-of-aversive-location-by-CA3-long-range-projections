% Figure 1

load('BigFatCluster.mat')
load('cclustID.mat')

%% A, B
% FOV, Example components and movement on running belt I-D mouse
mouse = 1;
session = 1;
nums = 1:15;
gain = [2 2];
helpers.ComponentsBehavePlot(CAIM,mouse,session,nums,gain)
%% C,D
% FOV, Example components and movement on running belt D-D mouse
mouse = 5;
session = 1;
nums = 1:15;
gain = [6 2];
helpers.ComponentsBehavePlot(CAIM,mouse,session,nums,gain)
%% E, F, G, H
% Event frequency of axons during run and rest periodes
session = 1:3;
helpers.EventFrequencyPlot(CAIM,cclustID,session)

%% I, J, K, L
% Place field, place preference Vector & Vector histogram plots
session = [3];
helpers.PlacePreferencePlot(CAIM,cclustID,session,[1 0]);

%% M
% Place preference vector length plot
figure
edges = 0:.02:1;
session = [1 2 3];
[PVvent,PVdors] = helpers.PlacePreferencePlot(CAIM,cclustID,session,[0 0]);

PVlv = histcounts(abs(PVvent),edges,"Normalization","probability");
PVld = histcounts(abs(PVdors),edges,"Normalization","probability");
plot(edges(1:end-1)+.01,cumsum(PVlv),"Color",[1 0 0])
hold on
plot(edges(1:end-1)+.01,cumsum(PVld),'Color',[0 0 0])
legend(["I-D" "D-D"],"Location","Northwest")
legend("Boxoff")
ylabel("Probability")
xlabel("Vector Length")
ylim([0 1])
