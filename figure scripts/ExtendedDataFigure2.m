
load('BigFatCluster.mat')
load('cclustID.mat')

%% C
session = [1 2 3 4 7 8];
numice = 1:4;
helpers.BehaviorPlot(CAIM,cclustID,session,numice)
sgtitle ('Behavior properties I-D')
numice = 5:8;
helpers.BehaviorPlot(CAIM,cclustID,session,numice)
sgtitle ('Behavior properties D-D')

%% D,E,F
% Event frequency of axons during run and rest periodes
for session = 1:3
    helpers.EventFrequencyPlot(CAIM,cclustID,session)
end