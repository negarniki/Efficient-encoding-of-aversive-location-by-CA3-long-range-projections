load('BigFatCluster.mat')
load('cclustID.mat')

%% A
session = [3 4 7 8];
numice = 1:8;
helpers.PlaceCoderPlot(CAIM,cclustID,session,numice);
%% B,C

figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[15 5 [15 7]],...
        'PaperUnits','centimeters',...
        'PaperSize', [15 7])

edges = 0:.02:1;
session = [3 4 7 8];
experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'Tn-1' 'Tn' 'P1'};
Style = {'-' '--' '-o' '-*'};
for i = 1:length(session)
    [PVvent,PVdors] = helpers.PlacePreferencePlot(CAIM,cclustID,session(i),[0 0]);
    PVlv = histcounts(abs(PVvent),edges,"Normalization","probability");
    PVld = histcounts(abs(PVdors),edges,"Normalization","probability");
    subplot(1,2,1)
    hold on
    plot(edges(1:end-1)+.01,cumsum(PVlv),Style{i},"Color",[1 0 0])
    subplot(1,2,2)
    hold on
    plot(edges(1:end-1)+.01,cumsum(PVld),Style{i},'Color',[0 0 0])
end

for i = 1:2
    subplot(1,2,i)
    xlim([0 1]);
    ylim([0 1]);
    legend(experiment(session),"Location","Northwest")
    legend("Boxoff")
    ylabel("Probability")
    xlabel("Vector Length")
end
%% D,E
session = [3 4 7 8];
helpers.PlacePreferencePlot(CAIM,cclustID,session,[0 1]);

%% F
% Place field density to running speed correlation
numexp = [3 4 7 8];
DoPlot = [1 0];
helpers.PlaceFieldDensityPlot(CAIM,cclustID,numexp,DoPlot);

%% G
% Place field density to running speed correlation
numexp = [1 2 3 4 7 8];
DoPlot = [0 1];
helpers.PlaceFieldDensityPlot(CAIM,cclustID,numexp,DoPlot);