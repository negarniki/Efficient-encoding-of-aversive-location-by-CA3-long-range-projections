load('BigFatCluster.mat')
load('cclustID.mat')

%% A,
session = [1 7 8];
helpers.SpeedPlot(CAIM,session,[1 0])

%% B
% Mean speed plot, linear and polar
session = [3 4 7 8];
helpers.SpeedPlot(CAIM,session,[0 1])

%% C
session = [3 7 8];
helpers.SpeedMinPlot(CAIM,session)