
%% Pool preference vectors of tracked components
DoPlot = [0 2];
pathname = 'C:\Users\martipof\IEECR Dropbox\Martin Pofahl\Matlab\NegarData\B1-B2\';
load([pathname 'BigFatClusterCorr.mat'])
numses = [1 2];
helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);
helpers.TrackedFOVPlot(CAIMcorr,cclustID)

pathname = 'C:\Users\martipof\IEECR Dropbox\Martin Pofahl\Matlab\NegarData\B3-T1\';
load([pathname 'BigFatClusterCorr.mat'])
numses = [1 2];
PV(:,:,2) = helpers.TrackedFeaturePlot(CAIMcorr,cclustID,numses,DoPlot);
helpers.TrackedFOVPlot(CAIMcorr,cclustID)