function [PVvent,PVdors] = PlacePreferencePlot(CAIM,cclustID,session,DoPlot)

experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'Tn-1' 'Tn' 'P1'};
mouseID = {'M259','M261', 'M270','M272','M262','M263','M271','M278'};

if DoPlot(1) == 1
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[5 1 [45 size(session,2)*4]],...
        'PaperUnits','centimeters',...
        'PaperSize', [45 size(session,2)*4])
end

PFventRet = [];
PFdorsRet = [];
PVvent = [];
PVventGr = [];
PVventMean= zeros(length(session),1);
PVdors = [];
PVdorsGr = [];
PVdorsMean = zeros(length(session),1);
PlacePrefVect = zeros(6,2,length(session));
numPlace = zeros(length(session),2);
numbin = 1500;
Alpha = linspace(pi,-pi,numbin*1);
Z1 = exp(Alpha(500)*1j);
Z2 = exp(Alpha(1000)*1j);
Z3 = exp(Alpha(1500)*1j);

for i = 1:length(session)
    dors = [];
    vent = [];
    PFvent = [];
    PCvent = [];
    PFdors = [];
    PCdors = [];
    for j = 1:size(CAIM,2)
        % Pool place preference vectors
        cclust = CAIM(session(i),j).cclust;
        isplace = cclust(:,cclustID.plcfld)>0 & cclust(:,cclustID.plcfldp)<=.05;
        PC(i,j) = sum(isplace);
        plccenter = cclust(isplace,cclustID.plcfld);
        plccenter = round(plccenter);
        plccenter(plccenter==0) = 1;
        plcvctang = cclust(isplace,cclustID.plcvctang);
        plcvctang = -(plcvctang-pi);
        plcvct = cclust(isplace,cclustID.plcvct);

        d = [];
        for k = 1:size(plccenter,1)        
            d = [d;plcvct(k)*exp(plcvctang(k)*1j)];
        end
        
        PF(i,j) = sum(d)/length(d);
        
        if j<5
            vent = [vent; d];
        else
            dors = [dors; d];
        end
                
        % Pool placefield maps
        plcfield = CAIM(session(i),j).plcfield; 
        plcfield = plcfield(isplace,:);
        plccenter = cclust(:,cclustID.plcfld);
        for k = 1:size(plcfield,1)
            plcfield(k,:) = plcfield(k,:)/max(plcfield(k,:));
        end
        if j<5
            PCvent = [PCvent; plccenter(isplace)];
            PFvent = [PFvent; plcfield(:,1:150)];
        else
            PCdors = [PCdors; plccenter(isplace)];
            PFdors = [PFdors; plcfield(:,1:150)];
        end
    end
    %%
    d = sum(vent)/length(vent); 
    PVvent = [PVvent; vent];   
    PVventGr(end+1:end+length(vent)) = session(i);%experiment(session(i));
    PVventMean(i) = abs(d);

    d = sum(dors)/length(dors); 
    PVdors = [PVdors; dors];   
    PVdorsGr(end+1:end+length(dors)) = session(i);%experiment(session(i));
    PVdorsMean(i) = abs(d);
    if DoPlot(1) == 1
        % Place field plot I-D projection
        subplot(length(session),6,(i-1)*6+1)
        [~,h] = sort(PCvent);
        PFvent = PFvent(h,:);
        imagesc(PFvent)
        colormap(hot)
        ax = gca;
    %     ax.YLim = [-.1 .2];
        ax.XLim = [0 150];
        
        temp = ax.Position;
        title(['I-D, Session: ' experiment{session(i)}])
        colorbar
        ax.Position = temp;
    
        % Polar plot of place preference vectors I-D projection
        subplot(length(session),6,(i-1)*6+2)
        
        for j = 1:length(vent)
            polarplot([0 real(-1j*log(vent(j)))],[0 abs(vent(j))],'color',[.5 .5 .5]);
            hold on
        end
        hold on
        
        polarplot([0 real(-1j*log(Z1))],[0 abs(Z1)],'color',[.1 .1 .8]);
        hold on
        polarplot([0 real(-1j*log(Z2))],[0 abs(Z2)],'color',[.1 .1 .8]);
        polarplot([0 real(-1j*log(Z3))],[0 abs(Z3)],'color',[.1 .1 .8]);
    
        grid off
        ax = gca;
        ax.RTickLabel = {};
    
        if CAIM(session(i),1).AP'*CAIM(session(i),1).stripe1.stimon>0
            DZ = exp(Alpha(500)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 60 120 240];
            ax.ThetaTickLabel = {'Zone2' 'AP' 'Zone1' 'Zone3' };
        elseif CAIM(session(i),1).AP'*CAIM(session(i),1).stripe2.stimon>0
            DZ = exp(Alpha(1000)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 120  240 320];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' 'AP'};
        elseif CAIM(session(i),1).AP'*CAIM(session(i),1).stripe3.stimon>0
            DZ = exp(Alpha(1500)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 120 180 320];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' 'AP' };
        else
            ax.ThetaTick = [0 120  240];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' };
        end
    
        % Mean place prefence vector plot I-D projection
        subplot(length(session),6,(i-1)*6+3)
        d = sum(vent)/length(vent);
        polarplot([0 real(-1j*log(d))],[0 abs(d)],'color',[0 0 0]);
        
        % Polar Histogram I-D projection
        polarhistogram(angle(vent),15,"Normalization","probability")
        hold on
        
        polarplot([0 real(-1j*log(Z1))],[0 abs(Z1)],'color',[.1 .1 .8]);
        hold on
        polarplot([0 real(-1j*log(Z2))],[0 abs(Z2)],'color',[.1 .1 .8]);
        polarplot([0 real(-1j*log(Z3))],[0 abs(Z3)],'color',[.1 .1 .8]);
    
        ax = gca;
        ax.RLim = [0 0.3];
        if CAIM(session(i),1).AP'*CAIM(session(i),1).stripe1.stimon>0
            DZ = exp(Alpha(500)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 60 120 240];
            ax.ThetaTickLabel = {'Zone2' 'AP' 'Zone1' 'Zone3' };
        elseif CAIM(session(i),1).AP'*CAIM(session(i),1).stripe2.stimon>0
            DZ = exp(Alpha(1000)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 120  240 320];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' 'AP'};
        elseif CAIM(session(i),1).AP'*CAIM(session(i),1).stripe3.stimon>0
            DZ = exp(Alpha(1500)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 120 180 320];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' 'AP' };
        else
            ax.ThetaTick = [0 120  240];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' };
        end
        
    %     p1 = hrtest(vent);
    %     [p1,T1] = hrtest(pi+angle(vent),10000);   
    %     p1 = round(p1,3);
    %     [p2,T2] = circ_rtest(pi+angle(vent));
    %     p2 = round(p2,3);
    %     disp(['p_{HR} = ' num2str(p1) ', p_{Ra} = ' num2str(p2)])
        %%
        % Place field plot D-D projection
        subplot(length(session),6,(i-1)*6+4)
        [~,h] = sort(PCdors);
        PFdors = PFdors(h,:);
        imagesc(PFdors)
        colormap(hot)
        ax = gca;
    %     ax.YLim = [-.1 .2];
        ax.XLim = [0 150];
        
        temp = ax.Position;
        title(['D-D'])
        colorbar
        ax.Position = temp;
    
        % Polar plot of place preference vectors D-D projection
        subplot(length(session),6,(i-1)*6+5)
        
        for j = 1:length(dors)
            polarplot([0 real(-1j*log(dors(j)))],[0 abs(dors(j))],'color',[.5 .5 .5]);
            hold on
        end
        hold on
        
        polarplot([0 real(-1j*log(Z1))],[0 abs(Z1)],'color',[.1 .1 .8]);
        hold on
        polarplot([0 real(-1j*log(Z2))],[0 abs(Z2)],'color',[.1 .1 .8]);
        polarplot([0 real(-1j*log(Z3))],[0 abs(Z3)],'color',[.1 .1 .8]);
    
        grid off
        ax = gca;
        ax.RTickLabel = {};
    
        if CAIM(session(i),5).AP'*CAIM(session(i),5).stripe1.stimon>0
            DZ = exp(Alpha(500)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 60 120 240];
            ax.ThetaTickLabel = {'Zone2' 'AP' 'Zone1' 'Zone3' };
        elseif CAIM(session(i),5).AP'*CAIM(session(i),5).stripe2.stimon>0
            DZ = exp(Alpha(1000)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 120  240 320];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' 'AP'};
        elseif CAIM(session(i),5).AP'*CAIM(session(i),5).stripe3.stimon>0
            DZ = exp(Alpha(1500)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 120 180 320];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' 'AP' };
        else
            ax.ThetaTick = [0 120  240];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' };
        end
        
        % Mean Polar plot D-D projection
        
        d = sum(dors)/length(dors);    
        polarplot([0 real(-1j*log(d))],[0 abs(d)],'color',[0 0 0]);
    
        
        % Polar histrogram of place preference vectors D-D projection
        subplot(length(session),6,(i-1)*6+6)
        
        polarhistogram(angle(dors),15,"Normalization","probability")
        hold on
        
        polarplot([0 real(-1j*log(Z1))],[0 abs(Z1)],'color',[.1 .1 .8]);
        hold on
        polarplot([0 real(-1j*log(Z2))],[0 abs(Z2)],'color',[.1 .1 .8]);
        polarplot([0 real(-1j*log(Z3))],[0 abs(Z3)],'color',[.1 .1 .8]);
    
        ax = gca;
        ax.RLim = [0 0.3];
    
        if CAIM(session(i),5).AP'*CAIM(session(i),5).stripe1.stimon>0
            DZ = exp(Alpha(500)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 60 120 240];
            ax.ThetaTickLabel = {'Zone2' 'AP' 'Zone1' 'Zone3' };
        elseif CAIM(session(i),5).AP'*CAIM(session(i),5).stripe2.stimon>0
            DZ = exp(Alpha(1000)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 120  240 320];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' 'AP'};
        elseif CAIM(session(i),5).AP'*CAIM(session(i),5).stripe3.stimon>0
            DZ = exp(Alpha(1500)*1j);
            polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
            ax.ThetaTick = [0 120 180 320];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' 'AP' };
        else
            ax.ThetaTick = [0 120  240];
            ax.ThetaTickLabel = {'Zone2'  'Zone1' 'Zone3' };
        end
    
    %     [p1,T1] = hrtest(pi+angle(dors));
    %     p1 = round(p1,3);
    %     [p2,T2]= circ_rtest(pi+angle(dors));
    %     p2 = round(p2,3);
    %     title(['p_{HR} = ' num2str(p1) ', p_{Ra} = ' num2str(p2)]);
    end
    %% Pool vectors for overview plot
    PlacePrefVect(1,:,i) = [abs(mean(vent)) abs(mean(dors))];
    % Std of length
    PlacePrefVect(2,:,i) = [std(abs(vent)) std(abs(dors))];
    % SEM of length
    PlacePrefVect(3,:,i) = PlacePrefVect(2,:,i)./sqrt([length(vent) length(dors)]);
    % Angle difference to AP
    dist = [(mean(vent)) (mean(dors))];
    shift = pi;%angle(DZ);%
    shift = exp(shift*1j);
    dist = dist./shift;
    dist = angle(dist);
    dist = -dist;
    dist = dist/pi*180;
    PlacePrefVect(4,:,i) = dist;
    % STD of Angles
    PlacePrefVect(5,:,i) = [std(vent) std(dors)]/pi*180;
    PlacePrefVect(6,:,i) = PlacePrefVect(5,:,i)./sqrt([length(vent) length(dors)]);
    %%
    numPlace(i,:) = [length(vent) length(dors)];
end

if DoPlot(2) == 1
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[15 5 [15 7]],...
        'PaperUnits','centimeters',...
        'PaperSize', [15 7])
    
    x = 1:size(PlacePrefVect,3);
    y1 = PlacePrefVect(1,1,:);
    y1 = permute(y1,[3 2 1]);
    y1err = PlacePrefVect(3,1,:);
    y1err = permute(y1err,[3 2 1]);
    y2 = PlacePrefVect(1,2,:);
    y2 = permute(y2,[3 2 1]);
    y2err = PlacePrefVect(3,2,:);
    y2err = permute(y2err,[3 2 1]);
    
    subplot(1,2,1)
    
    hold on
    errorbar(x,y1,y1err,'Linewidth',1,'Marker','.','MarkerSize',20,'Color','red')
    
    errorbar(x,y2,y2err,'Linewidth',1,'Marker','*','MarkerSize',10,'Color','black')
    % scatter(x,y2)
    
    ax = gca;
    ax.XTick = x;
    ax.XTickLabel = experiment(session);
    ax.XLim = [x(1)-1 x(end)+1];
    ylabel('Average vector length')
    
    y1 = PlacePrefVect(4,1,:);
    y1 = permute(y1,[3 2 1]);
    y1err = PlacePrefVect(6,1,:);
    y1err = permute(y1err,[3 2 1]);
    y2 = PlacePrefVect(4,2,:);
    y2 = permute(y2,[3 2 1]);
    y2err = PlacePrefVect(6,2,:);
    y2err = permute(y2err,[3 2 1]);
    
    subplot(1,2,2)
    
    hold on
    errorbar(x,y1,y1err,'Linewidth',1,'Marker','.','MarkerSize',20,'Color','red')
    
    errorbar(x,y2,y2err,'Linewidth',1,'Marker','*','MarkerSize',10,'Color','black')
    % scatter(x,y2)
    
    ax = gca;
    ax.XTick = x;
    ax.XTickLabel = experiment(session);
    ax.XLim = [x(1)-1 x(end)+1];
    ax.YLim = [-180 180];
    ylabel('Average angular distance')
    
    legend({'I-D' 'D-D'})
    legend('boxoff')
end
end