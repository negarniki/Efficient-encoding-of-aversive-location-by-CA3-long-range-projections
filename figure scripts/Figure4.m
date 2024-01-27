% Figure 4

% For violin plots, please download the following package:
% https://se.mathworks.com/matlabcentral/fileexchange/45134-violin-plot
% Otherwise data will be presened in BoxPlots

DoPlot = [1 0];
load(['B1-B2.mat'])
numses = [1 2];
[plcVentral1, plcDorsal1,DNrestActFract1,pooledf(1)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,0);

load(['B3-T1.mat'])
numses = [1 2];
[plcVentral2, plcDorsal2,DNrestActFract2,pooledf(2)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,0);

load([pathname 'B3-Tn.mat'])
numses = [1 2];
[plcVentral3, plcDorsal3,DNrestActFract3,pooledf(3)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,0);

load([pathname 'Tn-P.mat'])
numses = [1 2];
[plcVentral4, plcDorsal4,DNrestActFract4,pooledf(4)] = helpers.PlaceFieldShiftPlot(CAIMcorr,cclustID,numses,0);

%% C
ftsz = 8;
figSize=[21 7];
figure('color',[1 1 1],...
      'renderer','painters',...
      'visible','on',...
      'Units','centimeters',...
      'position',[20 5 figSize ],...
      'PaperUnits','centimeters',...
      'PaperSize', figSize )


ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
ax.XLim = [0 3];
ax.YLim = [0 1.5];
sgtitle ('Run-Rest activity - reference day i-D','fontsize', ftsz)



subplot(1,3,1)
diffLen2 =  length((pooledf(1).fVentral(:,2)))-length((pooledf(1).fDNrefVentral(:,2)));
diffLen3 =  length((pooledf(1).fVentral(:,2)))-length((pooledf(1).fDSrefVentral(:,2)));

% normalize and group the values for DNs and the other catagory
grouped1 = [
    [(pooledf(1).fDNrefVentral(:,2)-pooledf(1).fDNrefVentral(:,3))./pooledf(1).fDNrefVentral(:,1);nan(diffLen2, 1)],...
    [(pooledf(1).fDSrefVentral(:,2)-pooledf(1).fDSrefVentral(:,3))./pooledf(1).fDSrefVentral(:,1);nan(diffLen3, 1)]]; % group the values 

if exist("violin","file") == 2
    [h,L,MX,MED,bw] = violin (grouped1); % check the bandwidth
    bw = mean(bw); % make a mean of 3 bws
    violin (grouped1, bw,'facecolor', [1 1 0;.0 .0 .5;0 0.3 0.1]);
else
    boxplot(grouped1)
end

ax = gca;
ax.XTick = 1:3;
ax.XTickLabelRotation = 23;
ax.XTickLabel = {'DN' 'Others'};
ax.LineWidth = 1;
ax.FontSize = ftsz-2; 
ax.YLim = [-5 6];
ax.Box = "off";
plot([ax.XLim],[0 0],"--","color",[.5 .5 .5])
ylabel ('Firing Frequency Diff', 'fontsize',ftsz-2)
title("B_1-B_2 (B_1)")

subplot(1,3,2)
diffLen2 =  length((pooledf(3).fVentral(:,2)))-length((pooledf(3).fDNrefVentral(:,2)));
diffLen3 =  length((pooledf(3).fVentral(:,2)))-length((pooledf(3).fDSrefVentral(:,2)));

% normalize and group the values for DNs and the other catagory
grouped3 = [
    [(pooledf(3).fDNrefVentral(:,2)-pooledf(3).fDNrefVentral(:,3))./pooledf(3).fDNrefVentral(:,1);nan(diffLen2, 1)],...
    [(pooledf(3).fDSrefVentral(:,2)-pooledf(3).fDSrefVentral(:,3))./pooledf(3).fDSrefVentral(:,1);nan(diffLen3, 1)]];

if exist("violin","file") == 2
    [h,L,MX,MED,bw]= violin (grouped3); % check the bandwidth
    bw = mean(bw); % make a mean of 3 bws
    violin (grouped3, bw,'facecolor', [1 1 0;.0 .0 .5;0 0.3 0.1]); % plot again with the same bw
else
    boxplot(grouped3)
end

ax = gca;
ax.XTick = 1:3;
ax.XTickLabelRotation = 23;
ax.XTickLabel = {'DN' 'Others'};
ax.LineWidth = 1;
ax.FontSize = ftsz-2; 
ax.YLim = [-5 6];
ax.Box = "off";
plot([ax.XLim],[0 0],"--","color",[.5 .5 .5])
ylabel ('Firing Frequency Diff', 'fontsize',ftsz-2)
title("B_3-T_n (B_3)")

subplot(1,3,3)
% diffLen1 =  length((pooledf(4).fVentral(:,2)))-length((pooledf(4).fPCrefVentral(:,2)));
diffLen2 =  length((pooledf(4).fVentral(:,2)))-length((pooledf(4).fDNrefVentral(:,2)));
diffLen3 =  length((pooledf(4).fVentral(:,2)))-length((pooledf(4).fDSrefVentral(:,2)));

% normalize and group the values for DNs and the other catagory
grouped4 = [
    [(pooledf(4).fDNrefVentral(:,2)-pooledf(4).fDNrefVentral(:,3))./pooledf(4).fDNrefVentral(:,1);nan(diffLen2, 1)],...
    [(pooledf(4).fDSrefVentral(:,2)-pooledf(4).fDSrefVentral(:,3))./pooledf(4).fDSrefVentral(:,1);nan(diffLen3, 1)]];

if exist("violin","file") == 2
    [h,L,MX,MED,bw] = violin (grouped4); % check the bandwidth
    bw = mean(bw); % make a mean of 3 bws
    violin (grouped4, bw,'facecolor', [1 1 0;.0 .0 .5;0 0.3 0.1]); % plot again with the same bw
else
    boxplot(grouped4)
end

ax = gca;
ax.XTick = 1:3;
ax.XTickLabelRotation = 23;
ax.XTickLabel = {'DN' 'Others'};
ax.LineWidth = 1;
ax.FontSize = ftsz-2; 
ax.YLim = [-5 6];
ax.Box = "off";
plot([ax.XLim],[0 0],"--","color",[.5 .5 .5])
ylabel ('Firing Frequency Diff', 'fontsize',ftsz-2)
title("T_n-Pr (T_n)")

%%
figure('color',[1 1 1],...
      'renderer','painters',...
      'visible','on',...
      'Units','centimeters',...
      'position',[20 5 figSize ],...
      'PaperUnits','centimeters',...
      'PaperSize', figSize )


ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
ax.XLim = [0 3];
ax.YLim = [0 1.5];
% ax.LineWidth = lnwd;
sgtitle ('Run-Rest activity - reference day D-D','fontsize', ftsz)

subplot(1,3,1)
% diffLen1 =  length((pooledf(1).fDorsal(:,2)))-length((pooledf(1).fPCrefDorsal(:,2))); % PCs
diffLen2 = length((pooledf(1).fDorsal(:,2)))-length((pooledf(1).fDNrefDorsal(:,2))); % DN
diffLen3 = length((pooledf(1).fDorsal(:,2)))-length((pooledf(1).fDSrefDorsal(:,2)));%other


% normalize and group the values for DNs and the other catagory 
grouped1 = [[(pooledf(1).fDNrefDorsal(:,2)-pooledf(1).fDNrefDorsal(:,3))./pooledf(1).fDNrefDorsal(:,1); nan(diffLen2, 1)],...
    [(pooledf(1).fDSrefDorsal(:,2)-pooledf(1).fDSrefDorsal(:,3))./pooledf(1).fDSrefDorsal(:,1);nan(diffLen3, 1)]]; 

if exist("violin","file") == 2
    [h,L,MX,MED,bw] = violin (grouped1); % check the bandwidth
    bw = mean(bw); % make a mean of 3 bws
    violin (grouped1, bw,'facecolor', [1 1 0;.0 .0 .5;0 0.3 0.1]);
else
    boxplot(grouped1)
end

ax = gca;
ax.XTick = 1:3;
ax.XTickLabelRotation = 23;
ax.XTickLabel = {'DN' 'Others'};
ax.LineWidth = 1;
ax.FontSize = ftsz-2; 
ax.YLim = [-6 12];
ax.Box = "off";
plot([ax.XLim],[0 0],"--","color",[.5 .5 .5])
ylabel ('Firing Frequency Diff', 'fontsize',ftsz-2)
title("B_1-B_2 (B_1)")

subplot(1,3,2)
% diffLen1 =  length((pooledf(3).fDorsal(:,2)))-length((pooledf(3).fPCrefDorsal(:,2)));
diffLen2 =  length((pooledf(3).fDorsal(:,2)))-length((pooledf(3).fDNrefDorsal(:,2)));
diffLen3 =  length((pooledf(3).fDorsal(:,2)))-length((pooledf(3).fDSrefDorsal(:,2)));

% Normalize and group the values for DNs and others
grouped3 = [
    [(pooledf(3).fDNrefDorsal(:,2)-pooledf(3).fDNrefDorsal(:,3))./pooledf(3).fDNrefDorsal(:,1);nan(diffLen2, 1)],...
    [(pooledf(3).fDSrefDorsal(:,2)-pooledf(3).fDSrefDorsal(:,3))./pooledf(3).fDSrefDorsal(:,1);nan(diffLen3, 1)]]; 

if exist("violin","file") == 2
    [h,L,MX,MED,bw]= violin (grouped3); % check the bandwidth
    bw = mean(bw); % make a mean of 3 bws
    violin (grouped3, bw,'facecolor', [1 1 0;.0 .0 .5;0 0.3 0.1]); % plot again with the same bw
else
    boxplot(grouped3)
end

ax = gca;
ax.XTick = 1:3;
ax.XTickLabelRotation = 23;
ax.XTickLabel = {'DN' 'Others'};
ax.LineWidth = 1;
ax.FontSize = ftsz-2; 
ax.YLim = [-6 12];
ax.Box = "off";
plot([ax.XLim],[0 0],"--","color",[.5 .5 .5])
ylabel ('Firing Frequency Diff', 'fontsize',ftsz-2)
title("B_3-T_n (B_3)")


subplot(1,3,3)
diffLen2 =  length((pooledf(4).fDorsal(:,2)))-length((pooledf(4).fDNrefDorsal(:,2)));
diffLen3 =  length((pooledf(4).fDorsal(:,2)))-length((pooledf(4).fDSrefDorsal(:,2)));

% Normalize and group the values for DNs and others
grouped4 = [
    [(pooledf(4).fDNrefDorsal(:,2)-pooledf(4).fDNrefDorsal(:,3))./pooledf(4).fDNrefDorsal(:,1);nan(diffLen2, 1)],...
    [(pooledf(4).fDSrefDorsal(:,2)-pooledf(4).fDSrefDorsal(:,3))./pooledf(4).fDSrefDorsal(:,1);nan(diffLen3, 1)]]; 

if exist("violin","file") == 2
    [h,L,MX,MED,bw] = violin (grouped4); % check the bandwidth
    bw = mean(bw); % make a mean of 3 bws
    violin (grouped4, bw,'facecolor', [1 1 0;.0 .0 .5;0 0.3 0.1]); % plot again with the same bw
else
    boxplot(grouped4)
end

ax = gca;
ax.XTick = 1:3;
ax.XTickLabelRotation = 23;
ax.XTickLabel = {'DN' 'Others'};
ax.LineWidth = 1;
ax.FontSize = ftsz-2; 
ax.YLim = [-6 12];
ax.Box = "off";
plot([ax.XLim],[0 0],"--","color",[.5 .5 .5])
ylabel ('Firing Frequency Diff', 'fontsize',ftsz-2)
title("T_n-Pr (T_n)")