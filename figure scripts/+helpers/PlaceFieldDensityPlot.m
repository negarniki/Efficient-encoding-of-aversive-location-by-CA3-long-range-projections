function PlaceFieldDensityPlot(CAIM,cclustID,numexp,DoPlot)
%% Speed correlation using binned spatial data
% This function plots speed values and placefield density averaged accros
% all mice within each group for the experimental sessiond specified.
% It calculates pearsons' correlation coefficients and creates a grouped
% boxplot of these.
% For the grouped boxplots please download: Adam Danz (2024). boxplotGroup (https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup), MATLAB Central File Exchange. Retrieved January 22, 2024.


experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'TN-1' 'TN' 'P1' };

if DoPlot(1) == 1
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[3 3 [15 4*length(numexp)]],...
        'PaperUnits','centimeters',...
        'PaperSize', [15 4*length(numexp)])
end

x = 2.5:150/30:150;
PFVENT = zeros(length(numexp),30,size(CAIM,2)/2+1,3);
PFDORS = zeros(length(numexp),30,size(CAIM,2)/2+1,3);
SPEEDV = zeros(length(numexp),30,size(CAIM,2)/2+1,2);
SPEEDD = zeros(length(numexp),30,size(CAIM,2)/2+1,2);
PSr = zeros(length(numexp),size(CAIM,2));
PSrp = zeros(length(numexp),size(CAIM,2));
PSmeanr = zeros(length(numexp),2);
PSmeanrp = zeros(length(numexp),2);
YLimSpeed = [0.5 1.5];%[0 25]

for i = 1:length(numexp) 
    PFvent = [];
    PFdors = [];
    speedbinD = [];
    speedbinV = [];
    for j = 1:size(CAIM,2)
        
        cclust = CAIM(numexp(i),j).cclust;
        
        speedbin = CAIM(numexp(i),j).behave.speedbin;
        speedbin = speedbin(:,1:30);
        speedbin = speedbin*100;
        speedbin = speedbin/(100*mean(CAIM(numexp(i),j).behave.speedbin(2:end-1,1:30),'all'));
        for ii = 1:size(speedbin,1)
            speedbin(ii,:) = smooth(speedbin(ii,:),3);
        end
        
        if j<5
            speedbinV = [speedbinV;  speedbin];
        else
            speedbinD = [speedbinD;  speedbin];
        end
        
        speedbin = [mean(speedbin(2:end-1,:),1);
                std(speedbin(2:end-1,:),1)/sqrt(size(speedbin,1)-2)];
               
        isplace = cclust(:,cclustID.plcfld)>0 & cclust(:,cclustID.plcfldp)<=.05;
        plcvctang = cclust(isplace,cclustID.plcvctang);
        plcvctang(isnan(plcvctang)) = [];
        plcvctang(plcvctang<0) = plcvctang(plcvctang<0)+2*pi;
        bin = (2*pi)/size(speedbin,2);
        plccenter = zeros(1,size(speedbin,2));
        for ii = 1: size(speedbin,2)
            plccenter(ii) = sum(plcvctang>=(ii-1)*bin & plcvctang<ii*bin);
        end
%         temp = smooth(temp,3)';
        
        plcfield = CAIM(numexp(i),j).plcfield; 
        plcfield = plcfield(:,1:150);
        plcfield = plcfield(isplace,:);
%         plccenter = cclust(:,cclustID.plcfld);
%         plccenter = histcounts(plccenter(isplace)/10,0:5:150);

        
        bin = size(plcfield,2)/size(speedbin,2);
        temp = zeros(size(plcfield,1),size(speedbin,2));
        for k = 1:size(plcfield,1)
            for ii = 1: size(speedbin,2)
                temp(k,ii) = nanmean(plcfield(k,(ii-1)*bin+1:ii*bin));
            end
        end

        plcfield = temp;
        for k = 1:size(plcfield,1)
            plcfield(k,:) = plcfield(k,:)/max(plcfield(k,:));
        end
              
        if j<5
            PFvent = [PFvent; plcfield];
            PFVENT(i,:,j,1) = nanmean(plcfield,1);
            PFVENT(i,:,j,2) = nanstd(plcfield,1)/sqrt(size(plcfield,1));
            PFVENT(i,:,j,3) = plccenter;
            SPEEDV(i,:,j,1) = speedbin(1,:);
            SPEEDV(i,:,j,2) = speedbin(2,:);
        else
            PFdors = [PFdors; plcfield];
            PFDORS(i,:,j-4,1) = nanmean(plcfield,1);
            PFDORS(i,:,j-4,2) = nanstd(plcfield,1)/sqrt(size(plcfield,1));
            PFDORS(i,:,j-4,3) = plccenter;
            SPEEDD(i,:,j-4,1) = speedbin(1,:);
            SPEEDD(i,:,j-4,2) = speedbin(2,:);
        end
        
        plcfield = [nanmean(plcfield,1);
            nanstd(plcfield,[],1)/sqrt(size(plcfield,1))];

        [PSr(i,j),PSrp(i,j)] = corr(plcfield(1,:)',speedbin(1,:)');
        
    end
    
    [PSmeanr(i,1),PSmeanrp(i,1)] = corr(nanmean(PFvent,1)',nanmean(speedbinV)');
    [PSmeanr(i,2),PSmeanrp(i,2)] = corr(nanmean(PFdors,1)',nanmean(speedbinD)');
    
    PFVENT(i,:,5,1) = nanmean(PFvent,1);
    PFDORS(i,:,5,1) = nanmean(PFdors,1);
    SPEEDV(i,:,5,1) = nanmean(speedbinV,1);
    SPEEDD(i,:,5,1) = nanmean(speedbinD,1);
    
    PFVENT(i,:,5,2) = nanstd(PFvent)./sqrt(size(PFvent,1));
    PFDORS(i,:,5,2) = nanstd(PFdors)./sqrt(size(PFvent,1));
    SPEEDV(i,:,5,2) = nanstd(speedbinV)./sqrt(size(PFvent,1));
    SPEEDD(i,:,5,2) = nanstd(speedbinD)./sqrt(size(PFvent,1));
    
    if DoPlot(1) == 1
        subplot(length(numexp),2,(i-1)*2+1)
        y = nanmean(PFvent);
        yerr = nanstd(PFvent)/sqrt(size(PFvent,1));
        fill([x fliplr(x)],[y-yerr fliplr(y+yerr)],[.7 0 0],...
                'EdgeColor',[.7 0 0],...
                'EdgeAlpha',0,...
                'FaceAlpha',.5)
        hold on
        plot(x,y,'color',[.7 0 0],'linewidth',2,'LineStyle','-')
        title(['Session ' experiment{numexp(i)} ', I-D, r= ' num2str(round(PSmeanr(i,1),2))])
        ylim([0 .6])
        ylabel("PF density")
        xlabel('position on belt (cm)')
    
        yyaxis right
        y = nanmean(speedbinV);
        yerr = nanstd(speedbinV)/sqrt(size(speedbinV,1));
        fill([x fliplr(x)],[y-yerr fliplr(y+yerr)],[0 0 .7],...
                'EdgeColor',[0 0 .7],...
                'EdgeAlpha',0,...
                'FaceAlpha',.5)
        hold on
        plot(x,y,'color',[0 0 .7],'linewidth',2,'LineStyle','-')
        ylim(YLimSpeed)
        ax = gca;
        ax.YAxis(1).Color = [.7 0 0];
        ax.YAxis(2).Color = [0 0 .7];
        ylabel("normalized speed")
    
        subplot(length(numexp),2,(i-1)*2+2)
        y = nanmean(PFdors);
        yerr = nanstd(PFdors)/sqrt(size(PFdors,1));
        fill([x fliplr(x)],[y-yerr fliplr(y+yerr)],[.7 0 0],...
                'EdgeColor',[.7 0 0],...
                'EdgeAlpha',0,...
                'FaceAlpha',.5)
        hold on
        plot(x,y,'color',[.7 0 0],'linewidth',2,'LineStyle','-')
        title(['D-D, r= ' num2str(round(PSmeanr(i,2),2))])
        ylim([0 .6])
        ylabel("PF density")
        xlabel('position on belt (cm)')
    
        yyaxis right
        y = nanmean(speedbinD);
        yerr = nanstd(speedbinD)/sqrt(size(speedbinD,1));
        fill([x fliplr(x)],[y-yerr fliplr(y+yerr)],[0 0 .7],...
                'EdgeColor',[0 0 .7],...
                'EdgeAlpha',0,...
                'FaceAlpha',.5)
        hold on
        plot(x,y,'color',[0 0 .7],'linewidth',2,'LineStyle','-')
        ylim(YLimSpeed)
        ax = gca;
        ax.YAxis(1).Color = [.7 0 0];
        ax.YAxis(2).Color = [0 0 .7];
        ylabel("normalized speed")
    end
end

%%
if DoPlot(2) == 1
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[3 3 [20 10]],...
        'PaperUnits','centimeters',...
        'PaperSize', [20 10])
    xend = length(numexp)*3-2;
    x = [1:3:xend;1:3:xend;1:3:xend;1:3:xend];
    y = PSr(:,1:4)';
    scatter(x(:),y(:),'filled','MarkerFaceColor',[.5 .5 .5])
    hold on
    x = [2:3:xend+1;2:3:xend+1;2:3:xend+1;2:3:xend+1];
    y = PSr(:,5:8)';
    scatter(x(:),y(:),'filled','MarkerFaceColor',[.5 .5 .5])
    boxplotGroup({PSr(:,1:4)' PSr(:,5:8)'},'primaryLabels', {'ID' 'DD'},'secondaryLabels',experiment(numexp),'Color',[0 0 0]);
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    x = [1:3:xend];
    y = PSmeanr(:,1);
    scatter(x(:),y(:),'+','red')
    x = [2:3:xend+1];
    y = PSmeanr(:,2);
    scatter(x(:),y(:),'+','red')
    ylabel('r-coeff')
    title('ID vs DD speed-PFdensity-correlation')
    box off
end