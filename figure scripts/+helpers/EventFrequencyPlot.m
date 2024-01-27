function EventFrequencyPlot(CAIM,cclustID,session)

% plot Ca event frquencies

experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'Tn-1' 'Tn' 'P1'};
FrunVentral = [];
FrestVentral = [];
FrunDorsal = [];
FrestDorsal = []; 
PC = [];
PcV = [];
PcD = [];
for i = 1:length(session) 
    for k = 1:size(CAIM,2) 
        cclust = CAIM(session(i),k).cclust;       
        if k<5
            % access mean f for run
            FrunVentral = [FrunVentral; cclust(:,cclustID.meanf+1)];
            % access mean f for rest
            FrestVentral = [FrestVentral; cclust(:,cclustID.meanf+2)];
        else
%             access mean f for run
            FrunDorsal = [FrunDorsal; cclust(:,cclustID.meanf+1)];
            % access mean f for rest
            FrestDorsal = [FrestDorsal; cclust(:,cclustID.meanf+2)];
        end
        
        cclust = CAIM(session(i),k).cclust;
        % losonczy criterium
%         isplace = cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcvctp)<=.05;
        % Dombeck criterium
        isplace = cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcfldp)<=.05;
        % all cells 
%         isplace = cclust(:,cclustID.plcvct)>0%& cclust(:,cclustID.plcvctp)>=.05;
                     
        if k < 5
              PcV = [PcV;isplace];
        else
              PcD = [PcD;isplace];
        end    
    end   
end


%% BoxPlot mean_firing frequency 
ftsz = 16;
figSize=[15 15];
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
% ax.XTick = [1 2];
% ax.XTickLabel = {'rest' 'run'};
box off

plot(1:2,[FrunVentral  FrestVentral ]','color',[.5 .5 .5])

hold on
plot(3:4,[FrunDorsal FrestDorsal]','color',[.5 .5 .5])

diffLen =  length(FrunVentral)-length(FrunDorsal);
if diffLen>0
    grouped = [FrunVentral, FrestVentral, [FrunDorsal;nan(diffLen, 1)],  [FrestDorsal;nan(diffLen, 1)]];
else 
    grouped = [[FrunVentral; nan(diffLen, 1)], [FrestVentral;nan(diffLen, 1)], FrunDorsal,  FrestDorsal];
end

% figure
boxplot(grouped)

ax = gca;
ax.XTick = 1:4;
% ax.XTickLabel = {'run_ventral' 'run_Dorsal' 'rest_ventral' 'rest_Dorsal'};
xtickangle(23)
ax.XTickLabel = {'run ventral' 'rest ventral' 'run dorsal' 'rest dorsal'};
set(gca,'box', 'off', 'linewidth', 1,'fontsize',ftsz-2)%,'fontname','FreeSerif') 

ylabel ('Mean Firing Frequency', 'fontsize',ftsz-2)
title([experiment{session}], 'fontsize',ftsz)
box off

%% Violin plot

% For violin plots, please download the following package:
% https://se.mathworks.com/matlabcentral/fileexchange/45134-violin-plot
if exist("violin",'file') ==2
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
    % ax.XTick = [1 2];
    % ax.XTickLabel = {'rest' 'run'};
    box off
    plot(1:2,[FrunVentral  FrestVentral ]','color',[.5 .5 .5])
    
    hold on
    plot(3:4,[FrunDorsal FrestDorsal]','color',[.5 .5 .5])
    
    hold on
    
    
    violin(grouped, 'facecolor', [.4 .2 1])
    
    ax = gca;
    ax.XTick = 1:4;
    xtickangle(23)
    ax.XTickLabel = {'run ventral' 'rest ventral' 'run dorsal' 'rest dorsal'};
    set(gca,'box', 'off', 'linewidth', 1,'fontsize',ftsz-2)%,'fontname','FreeSerif') 
    
    ylabel ('Mean Firing Frequency', 'fontsize',ftsz-2)
    title([experiment{session}], 'fontsize',ftsz)
    box off
end


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
% ax.XTick = [1 2];
% ax.XTickLabel = {'rest' 'run'};
box off

% Plot the Ventral run and rest active fibers against each other and label the place coders

plot(FrestVentral,FrunVentral,'ok')
% plot(FrestDorsal,FrunDorsal,'ok')

xlabel('REST')
ylabel('RUN')
hold on
% figure (1)
plot(FrestVentral(FrunVentral>FrestVentral),FrunVentral(FrunVentral>FrestVentral),'or')
plot([0  max([FrunVentral;FrestVentral])],[0 max([FrunVentral;FrestVentral])],'--k')
% comment in for ventral 
plot(FrestVentral(find(PcV)),FrunVentral(find(PcV)),'b*')
xlim([0 12])
ylim([0 12])

title(['Mean Firing Frequency Ventral - ' [experiment{session}] ])
legend ('Rest','Run','Unity line','Place Coders');legend boxoff
% legend ('Rest','Run','Unity line');legend boxoff
axis square
box off 

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
% ax.XTick = [1 2];
% ax.XTickLabel = {'rest' 'run'};
box off

%Plot the Dorsal run and rest active fibers against each other and label the place coders
plot(FrestDorsal,FrunDorsal,'ok')

xlabel('REST')
ylabel('RUN')
hold on
% ([FrunDorsal;FrestDorsal]),'--k')

plot(FrestDorsal(FrunDorsal>FrestDorsal),FrunDorsal(FrunDorsal>FrestDorsal),'or')
plot([0  max([FrunDorsal;FrestDorsal])],[0 max([FrunDorsal;FrestDorsal])],'--k')
plot(FrestDorsal(find(PcD)),FrunDorsal(find(PcD)),'b*')

xlim([0 12])
ylim([0 12])

title (['Mean Firing Frequency Dorsal - ' [experiment{session}] ])
legend ('Rest','Run','Unity line','Place Coders');legend boxoff
% legend ('Rest','Run','Unity line');legend boxoff
axis square
box off
% print(gcf, '-dpdf', [pdfpathname 'Dors_Rest_Run.pdf'])


 

end