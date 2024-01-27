
%% A,B
% Plot the pplace preference Vectors of tracked feature coders for all transitions
DoPlot = 1;
load(['B1-B2.mat'])
numses = [1 2];
[plcVentral1, plcDorsal1,DNrestActFract1,pooledf(1)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,DoPlot);

load(['B3-T1.mat'])
numses = [1 2];
[plcVentral2, plcDorsal2,DNrestActFract2,pooledf(2)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,DoPlot);

load([pathname 'B3-Tn.mat'])
numses = [1 2];
[plcVentral3, plcDorsal3,DNrestActFract3,pooledf(3)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,DoPlot);

load([pathname 'Tn-P.mat'])
numses = [1 2];
[plcVentral4, plcDorsal4,DNrestActFract4,pooledf(4)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,DoPlot);

