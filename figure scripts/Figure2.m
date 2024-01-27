
load('BigFatCluster.mat')
load('cclustID.mat')

%% B, C
numexp = [3 4 7 8];
helpers.PlacePreferencePlot(CAIM,cclustID,numexp,[1 0]);
