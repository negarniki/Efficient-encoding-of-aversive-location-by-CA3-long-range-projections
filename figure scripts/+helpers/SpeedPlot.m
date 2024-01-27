function SpeedPlot(CAIM,session,DoPlot)
    
experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'TN-1' 'TN' 'P1' };
mouseID = {'M259','M261', 'M270','M272','M262','M263','M271','M278'};

SP = zeros(length(session),30);
SPstd = zeros(length(session),30);
SPmax = [];
SPref = [];
SPgr = cell(0);%[];
SPmouse = cell(0);
WT = zeros(length(session),30);
WTstd = zeros(length(session),30);
int = 9:10; % interval do check for speed differences
refint = 5:6;
refses = 3;
x = 2.5:5:150;
YLim = [5 35];

if DoPlot(1) == 1
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[5 1 [45 25]],...
        'PaperUnits','centimeters',...
        'PaperSize', [45 25])
end

for i = 1:length(session)
    waitbinD = [];
    waitbinV = [];
    speedbinD = [];
    speedbinV = [];
    for j = 1:size(CAIM,2)
        behave = CAIM(session(i),j).behave;
        
        if DoPlot(1) == 1
            subplot(length(session)*2,size(CAIM,2),2*(i-1)*(size(CAIM,2))+j)
            imagesc(x,size(behave.speedbin,2)-2,behave.speedbin(2:end-1,1:30))
            if j == 1
                title(['session ' experiment{session(i)} ', ' mouseID{j}])
            else
                title(mouseID{j})
            end
            
            subplot(length(session)*2,size(CAIM,2),2*(i-1)*(size(CAIM,2))+size(CAIM,2)+j)
            hold on
        end

        speedtemp = behave.speedbin(2:end-1,1:30);
        speedtemp = speedtemp*100;
%         speedtemp = speedtemp/(100*mean(CAIM(session(i),j).behave.speedbin(2:end-1,1:30),'all'));
        
        for ii = 1:size(speedtemp,1)
            temp = speedtemp(ii,:);
            temp = [temp(end-9:end) temp temp(1:10)];
            temp = smooth(temp,3);
%             temp = diff(temp);temp(end+1) = temp(end);
            speedtemp(ii,:) = temp(11:end-10);
        end
        
        if DoPlot(1) == 1
            y = nanmean(speedtemp,1);
            yy = nanstd(speedtemp,1)/sqrt(size(speedtemp,1));       
            fill([x fliplr(x)],[y-yy fliplr(y+yy)],[0 0 .8],...
            'EdgeColor',[0 0 .8],...
            'EdgeAlpha',0,...
            'FaceAlpha',.5)
            plot(x,y,'linewidth',1,'color',[0 0 .8])
            ylim(YLim)  
            xlabel(['position (cm)'])
            ylabel(['speed (cm/s)'])
        end

        if j<5
            speedbinV = [speedbinV;  speedtemp];
        else
            speedbinD = [speedbinD; speedtemp];
        end
        SPmouse(end+1:end+size(speedtemp,1),1) = mouseID(j); 
    end

    y = nanmean([speedbinV; speedbinD]);
    yy = nanstd([speedbinV;speedbinD])/sqrt(size([speedbinV;speedbinD],1));
    [~,minSpeed] = min(y(5:end-5)); 
    minSpeedD = (minSpeed+4)*5-2.5;
    meanSpeed = mean([speedbinV; speedbinD],'all');

    SP(i,:) = y;
    SPstd(i,:) = yy/sqrt(size(speedbinD,1)+size(speedbinV,1));
    WT(i,:) = nanmean([waitbinD; waitbinV]);
    WTstd(i,:) = nanstd([waitbinD; waitbinV])/sqrt(8);
    
   
%     if session(i) == 3 || session(i) == 5   || session(i) == 6 || session(i) == 10 || session(i) == 13 || session(i) == 14
%     SPmax = [SPmax; [[max(speedbinD(:,int),[],2); max(speedbinV(:,int),[],2)] [max(speedbinD(:,int+10),[],2); max(speedbinV(:,int+10),[],2)] [max(speedbinD(:,int+20),[],2); max(speedbinV(:,int+20),[],2)]]];
    SPmax = [SPmax; [[mean(speedbinV(:,int+1),2); mean(speedbinD(:,int+1),2)] [mean(speedbinV(:,int+10),2); mean(speedbinD(:,int+10),2)] [mean(speedbinV(:,int+20),2); mean(speedbinD(:,int+20),2)]]];
    SPref = [SPref; [mean(speedbinV(:,refint),2); mean(speedbinD(:,refint),2)]];
    SPgr(end+1:end+size(speedbinV,1)+size(speedbinD,1),1) = experiment(session(i));
%     end
end
 
%%
if DoPlot(2) == 1
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[20 5 [12 20]],...
        'PaperUnits','centimeters',...
        'PaperSize', [12 20])
    
    subplot(2,1,1)
    
    % put in what ever session you desire (These Sessions need to be loaded in
    % the sectionb above)
    
    num = 1:length(session);
    numcol = hsv(length(num));
    % numcol = zeros(length(num),3);
    p = [];
    x = 0:5:150;
    for i = 1:length(num)
        y = [SP(num(i),:) SP(num(i),1)];%+(length(num)-i)*8;
        yerr = [SPstd(num(i),:) SPstd(num(i),1)];
        fill([x fliplr(x)],[y-yerr fliplr(y+yerr)],numcol((i),:),...
                'EdgeColor',[0 .8 0],...
                'EdgeAlpha',0,...
                'FaceAlpha',.5)
            hold on
        p(i) = plot(x,y,'color',numcol((i),:));
    
        xlabel('position (cm)')
        ylabel('mean speed')
    end
    
    legend(p,experiment(session(num)))
    legend('boxoff')
    
    subplot(2,1,2)
    
    % num = [3 5 6];
    Alpha = linspace(pi,-pi,size(SP,2)+1);
    Z1 = exp(Alpha(11)*1j);
    Z2 = exp(Alpha(21)*1j);
    Z3 = exp(Alpha(31)*1j);
    
    RLim = [10 22];
    polarplot([0 real(-1j*log(Z1))],RLim,'color',[.1 .1 .8]);
    hold on
    polarplot([0 real(-1j*log(Z2))],RLim,'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(Z3))],RLim,'color',[.1 .1 .8]);
    
    for i = 1:length(num)
        polarplot(Alpha,[SP(num(i),:) SP(num(i),1)],'color',numcol((i),:))
        hold on
    end
    ax = gca;
    ax.ThetaTick = [0 120  240];
    ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' };
    ax.RLim = RLim;
end
end