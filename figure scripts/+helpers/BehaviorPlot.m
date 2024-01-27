function BehaviorPlot(CAIM,cclustID,session,numice)
% this code plots the basic behavioral properties 
    %%

    lnwd = 1;
    ftsz = 8;
    pupilsize = nan(3,length(session),length(numice));
    totdist = nan(length(session),length(numice));
    meanspeed = nan(length(session),length(numice));
    runtime = nan(length(session),length(numice));
    for i = 1:length(session)
        for j = 1:length(numice)           
            pupilsize(:,i,j) = CAIM(session(i),numice(j)).behave.pupilsize(1,:);
            pupilsize(:,i,j) = pupilsize(:,i,numice(j))./CAIM(session(i),numice(j)).behave.pupilsize(1,1);
            totdist(i,j) = CAIM(session(i),numice(j)).behave.totdist;
            meanspeed(i,j) = CAIM(session(i),numice(j)).behave.meanspeed;
            runtime(i,j) = sum(CAIM(session(i),numice(j)).behave.running==1)/length(CAIM(session(i),numice(j)).behave.running);
        end
    end
    
    %%
    figsize = [21 7];
    %plot total dist
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[20 5 figsize],...
        'PaperUnits','centimeters',...
        'PaperSize', figsize)
    
    subplot(1,3,1)
    
    x = [1 2 3 4]; 
    % put here again the right sessions
    y = totdist(1:3,:)/1000;
    y = y(:)';
    yy = totdist(4,:)/1000;
    yy = yy(:)';
    yyy = totdist(5,:)/1000;
    yyy = yyy(:)';
    yyyz = totdist(6,:)/1000;
    
    
    z = meanspeed(1:3,:);
    z = z(:)';
    zz = meanspeed(4,:);
    zz = zz(:)';
    zzz = meanspeed (5,:);
    zzz = zzz(:)'; 
    zzzz = meanspeed(6,:);
    % yyyz = yyyz(:)';
    
    
    yerr = [nanstd(y,0,2)./sqrt(size(y,2)); nanstd(yy,0,2)./sqrt(size(yy,2));nanstd(yyy,0,2)./sqrt(size(yyy,2));nanstd(yyyz,0,2)./sqrt(size(yyyz,2))];
    y = [nanmean(y,2); nanmean(yy,2);nanmean(yyy,2); nanmean(yyyz,2)];
    
    b = bar(x, y);
    
    hold on
    % b.FaceColor = 'flat';
    % b.LineStyle = 'none';
    % b.CData(1,:) = [.2 .2 .2];
    % b.CData(2,:) = [.5 .5 .5];
    
    errorbar(x,y,nan(size(yerr)),yerr,'.',...
                'Marker','none',...
                'Color',[0 0 0],...
                'LineWidth',lnwd)
            
    % title(experiment{5},'FontWeight','normal')
            
    ax = gca;
    ax.FontSize = ftsz-2;
    ax.YColor = [0 0 0];
    ax.XColor = [0 0 0];
    ax.Color = [1 1 1];
    ax.XLim = [0 5];
    ax.YLim = [0 max(y)+max(yerr)+1];
    ax.YLim = [0 50];
    ax.LineWidth = lnwd;
    ax.XTick = x;
    ax.XTickLabelRotation = -30;
    ax.XTickLabel = {'B 1-3' 'T1' 'Tn' 'P1'};
    ylabel('total distance (m)','FontSize',ftsz-2)
    box('off')
    hold off
    
    subplot(1,3,2)
    
    x = [1 2 3 4]; 
    
    y = meanspeed(1:3,:);
    y = y(:)';
    yy = meanspeed(4,:);
    yy = yy(:)';
    yyy = meanspeed (5,:);
    yyy = yyy(:)'; 
    yyyz = meanspeed(6,:);
    % yyyz = yyyz(:)';
    
    yerr = [nanstd(y,0,2)./sqrt(size(y,2)); nanstd(yy,0,2)./sqrt(size(yy,2)); nanstd(yyy,0,2)./sqrt(size(yyy,2));nanstd(yyyz,0,2)./sqrt(size(yyyz,2))];
    y = [nanmean(y,2); nanmean(yy,2); nanmean(yyy,2);nanmean(yyyz,2)];
    
    b = bar(x, y);
    
    hold on
    % b.FaceColor = 'flat';
    % b.LineStyle = 'none';
    % b.CData(1,:) = [.2 .2 .2];
    % b.CData(2,:) = [.5 .5 .5];
    
    errorbar(x,y,nan(size(yerr)),yerr,'.',...
                'Marker','none',...
                'Color',[0 0 0],...
                'LineWidth',lnwd)
            
    % title(experiment{5},'FontWeight','normal')
            
    ax = gca;
    ax.FontSize = ftsz-2;
    ax.YColor = [0 0 0];
    ax.XColor = [0 0 0];
    ax.Color = [1 1 1];
    ax.XLim = [0 5];
    % ax.YLim = [0 max(y)+max(yerr)+1];
    ax.YLim = [0 15];
    ax.LineWidth = lnwd;
    ax.XTick = x;
    ax.XTickLabelRotation = -30;
    % ax.XTickLabel = {'B 1-2' 'T 1-2' 'P 1-2' 'E 1-2'};
    ax.XTickLabel = {'B 1-3' 'T1' 'Tn' 'P1'};
    ylabel('mean speed (cm s^{-1})','FontSize',ftsz-2)
    box('off')
    hold off
    
    subplot(1,3,3)
    
    x = [1 2 3 4]; 
    
    y = runtime(1:3,:)*100;
    y = y(:)';
    yy = runtime(4,:)*100;
    yy = yy(:)';
    yyy = runtime(5,:)*100;
    yyy = yyy(:)';
    yyyz = runtime(6,:)*100;
    % yyy = yyy(:)';
    
    yerr = [nanstd(y,0,2)./sqrt(size(y,2)); nanstd(yy,0,2)./sqrt(size(yyy,2)); nanstd(yyy,0,2)./sqrt(size(yyy,2));nanstd(yyyz,0,2)./sqrt(size(yyyz,2)) ];
    y = [nanmean(y,2); nanmean(yy,2);nanmean(yyy,2);nanmean(yyyz,2)];
    
    b = bar(x, y);
    
    hold on
    % b.FaceColor = 'flat';
    % b.LineStyle = 'none';
    
    % b.CData(1,:) = [.2 .2 .2];
    % b.CData(2,:) = [.5 .5 .5];
    
    errorbar(x,y,nan(size(yerr)),yerr,'.',...
                'Marker','none',...
                'Color',[0 0 0],...
                'LineWidth',lnwd)
            
    % title(experiment{5},'FontWeight','normal')
            
    ax = gca;
    ax.FontSize = ftsz-2;
    ax.YColor = [0 0 0];
    ax.XColor = [0 0 0];
    ax.Color = [1 1 1];
    ax.XLim = [0 5];
    % ax.YLim = [0 max(y)+max(yerr)+1];
    ax.YLim = [0 40];
    ax.LineWidth = lnwd;
    ax.XTick = x;
    ax.XTickLabelRotation = -30;
    ax.XTickLabel = {'B 1-3' 'T1' 'Tn' 'P1'};
    % ax.XTickLabel = {'B 1-2' 'T 1-2' 'P 1-2' 'E 1-2'};
    ylabel('fraction of locomotion time (%)','FontSize',ftsz-2)
    box('off')
    hold off
    % print(gcf, '-dpdf', [pdfpathname 'runtime_i-d_old.pdf'])

end
   