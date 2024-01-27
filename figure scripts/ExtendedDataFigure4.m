
load('BigFatCluster.mat')
load('cclustID.mat')

%% A, B
numexp = [1 2 3];
helpers.PlacePreferencePlot(CAIM,cclustID,numexp,[1 0]);

%% C
numexp = [1 2 3];
helpers.PlaceFieldDensityPlot(CAIM,cclustID,numexp,[1 0]);
