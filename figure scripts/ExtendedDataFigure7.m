load('BigFatCluster.mat')
load('cclustID.mat')

%%
numexp = [3 4 7 8];
experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'TN-1' 'TN' 'P1' };
mouseID = {'M259','M261', 'M270','M272','M262','M263','M271','M278'};

SP = zeros(length(numexp),30);
SPstd = zeros(length(numexp),30);
SPmax = [];
SPref = [];
SPgr = cell(0);%[];
SPmouse = cell(0);
WT = zeros(length(numexp),30);
WTstd = zeros(length(numexp),30);
int = 9:10; % interval do check for speed differences
refint = 5:6;
refses = 3;
x = 2.5:5:150;
YLim = [5 35];

for i = 1:length(numexp)
    waitbinD = [];
    waitbinV = [];
    speedbinD = [];
    speedbinV = [];
    for j = 1:size(CAIM,2)
        behave = CAIM(numexp(i),j).behave;
        
        speedtemp = behave.speedbin(2:end-1,1:30);
        speedtemp = speedtemp*100;
        speedtemp = speedtemp/(100*mean(CAIM(numexp(i),j).behave.speedbin(2:end-1,1:30),'all'));
        waittemp = behave.waitbin(2:end-1,1:30);
        waittemp = mat2gray(waittemp);
        for ii = 1:size(speedtemp,1)
            temp = speedtemp(ii,:);
            temp = [temp(end-9:end) temp temp(1:10)];
            temp = smooth(temp,3);
            speedtemp(ii,:) = temp(11:end-10);
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
    
   
%     if numexp(i) == 3 || numexp(i) == 5   || numexp(i) == 6 || numexp(i) == 10 || numexp(i) == 13 || numexp(i) == 14
%     SPmax = [SPmax; [[max(speedbinD(:,int),[],2); max(speedbinV(:,int),[],2)] [max(speedbinD(:,int+10),[],2); max(speedbinV(:,int+10),[],2)] [max(speedbinD(:,int+20),[],2); max(speedbinV(:,int+20),[],2)]]];
    SPmax = [SPmax; [[mean(speedbinV(:,int+1),2); mean(speedbinD(:,int+1),2)] [mean(speedbinV(:,int+10),2); mean(speedbinD(:,int+10),2)] [mean(speedbinV(:,int+20),2); mean(speedbinD(:,int+20),2)]]];
    SPref = [SPref; [mean(speedbinV(:,refint),2); mean(speedbinD(:,refint),2)]];
    SPgr(end+1:end+size(speedbinV,1)+size(speedbinD,1),1) = experiment(numexp(i));
%     end
end

%% B, C, D, E
% grouped boxplots every individual mouse AP vs erf RW vs ref

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[20 5 [25 25]],...
    'PaperUnits','centimeters',...
    'PaperSize', [25 25])

stripe = [1 2 3];
strcomp = {'I-AP' 'II' 'III-RW' 'ref'};
YLim = [0 2.1];
grp = zeros(length(stripe),1);
grn = zeros(length(stripe),1);
ttp = zeros(length(stripe),length(mouseID));
ttn = zeros(length(stripe),length(mouseID));
j = 3;
for k = 1:length(stripe)
    subplot(4,5,(k-1)*5+(1:4))    
    x1 = [];
    x2 = [];
    for i = 1:length(mouseID)   
    
        temp = ismember(SPmouse,mouseID(i)) & ismember(SPgr,experiment(numexp(j)));
        x1(1:sum(temp),i) = SPmax(temp,stripe(k));
        x2(1:sum(temp),i) = SPref(temp,1);
        [~,ttp(k,i)] = ttest(x1(1:sum(temp),i),x2(1:sum(temp),i));
        ttn(k,i) = sum(temp);
        plot([-2 -1]+i*3,[x1(1:sum(temp),i) x2(1:sum(temp),i)],'color',[.8 .8 .8])
        hold on
        if ttp(k,i)<.05;text(-1.5+i*3,YLim(2)-.1,'*');end
    end
    x1(x1==0) = nan;
    x2(x2==0) = nan;
    x = [{x1} {x2}];
    boxplotGroup(x,'primaryLabels', [strcomp(k) strcomp(end)],'secondaryLabels',mouseID);
    title([experiment{numexp(j)} ', ' strcomp{k} ' vs. ' strcomp{end}] )
    box off
    hold off
    ylim(YLim)
    ylabel('running speed/mean')
    temppos = get(gca,'position');

    subplot(4,5,k*5)
    temp = ismember(SPgr,experiment(numexp(j)));
    x = [];
    x(1:sum(temp),1) = SPmax(temp,stripe(k));
    x(1:sum(temp),2) = SPref(temp,1);
    plot(1:2,x','color',[.8 .8 .8])
    hold on
    boxplot(x,{strcomp{k} strcomp{end}})
    box off
    temppos(2,:) = get(gca,'position');
    set(gca,'position',[temppos(2,1) temppos(1,2) temppos(2,3) temppos(1,4)])
    ylim(YLim)
    ylabel('running speed/mean')
    [grp(k),tbl,stats] = anova1(x,[],"off");
    grn(k) = size(x,1);
    title(['p = ' num2str(grp(k))])
end

YLim = [-1 1];
tbl = cell(1,4);
for j = 1:4
    subplot(4,4,12+j)
    temp = ismember(SPgr,experiment(numexp(j)));
    x = SPmax - SPref;
    x = x(temp,:);
    plot(1:3,x','color',[.8 .8 .8])
    hold on
    boxplot(x,{'I-AP' 'II' 'III-RW'})    
    ylabel('\Delta running speed')
    box off
    ylim(YLim)
    
    [p,tbl{j},stats] = anova1(x,{'I-AP' 'II' 'III-RW'},"off");
    [c,m,h,nms] = multcompare(stats,'Estimate','column','CType','bonferroni','Display','off');
    if p<.05
        if c(1,end)<.05
            plot([1.1 1.9], [YLim(2)-.3 YLim(2)-.3],'color',[0 0 0])
            text(1.5,YLim(2)-.2,'*')
        end
        if c(2,end)<.05
            plot([1.1 2.9], [YLim(2)-.1 YLim(2)-.1],'color',[0 0 0])
            text(2,YLim(2),'*');
        end
        if c(3,end)<.05
            plot([2.1 2.9], [YLim(2)-.3 YLim(2)-.3],'color',[0 0 0])
            text(2.5,YLim(2)-.2,'*');
        end
    end
    pmult(:,j)=c(:,end);
    title([experiment{numexp(j)} ', p = ' num2str(p)])
end

%% F,G
% Correlating changes of speed with changes of 'place fields/Placecode density ENTIRE ZONE

shift = -0;
int = 1:10;
refint = 11:20;
speedint = 9:10;
speedref = 12:14;
readout = 3;% 1 = PF density, 3 = place field centers
xrange = [-1 .5];%[-15 5]
yrange = [-2.5 2];%[-.2 .5];%

PFeff = zeros(length(numexp),4,4);
SPeff = zeros(length(numexp),4,4);
for i = 1:length(numexp) % Session
    for j = 1:4 % Mouse
        PFeff(i,j,1) = mean(PFVENT(i,int+shift,j,readout))-mean(PFVENT(i,refint,j,readout));
        SPeff(i,j,1) = mean(SPEEDV(i,speedint+1,j,1))-mean(SPEEDV(i,speedref,j,1));
        PFeff(i,j,2) = mean(PFVENT(i,int+20+shift,j,readout))-mean(PFVENT(i,refint,j,readout));
        SPeff(i,j,2) = mean(SPEEDV(i,speedint+20,j,1))-mean(SPEEDV(i,speedref,j,1));

        PFeff(i,j,3) = mean(PFDORS(i,int+shift,j,readout))-mean(PFDORS(i,refint,j,readout));
        SPeff(i,j,3) = mean(SPEEDD(i,speedint+1,j,1))-mean(SPEEDD(i,speedref,j,1));
        PFeff(i,j,4) = mean(PFDORS(i,int+20+shift,j,readout))-mean(PFDORS(i,refint,j,readout));
        SPeff(i,j,4) = mean(SPEEDD(i,speedint+20,j,1))-mean(SPEEDD(i,speedref,j,1));
    end
end
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[3 3 [25 15]],...
    'PaperUnits','centimeters',...
    'PaperSize', [25 15])

subplot(1,2,1)
refses = 1:4;
x = SPeff(refses,:,1);
y = PFeff(refses,:,1);
[r(1),p(1)] = corr(x(:),y(:));
scatter(x(:),y(:),'filled','MarkerFaceColor',[.8 0 0])
hold on
x = SPeff(refses,:,2);
y = PFeff(refses,:,2);
[r(2),p(2)] = corr(x(:),y(:));
scatter(x(:),y(:),'filled','MarkerFaceColor',[0 .8 0])
xlim(xrange)
ylim(yrange)
hold off
xlabel('\Delta Running Speed (cm/s)')
% ylabel('\Delta Place Field Density (AU)')
ylabel('\Delta PV vectors in bin (AU)')
title({['I-D, ' experiment{numexp(refses)} ', Response Shift: ' num2str(shift*5) ' cm']...
    ['AP: r = ' num2str(round(r(1),2)) ', p = ' num2str(round(p(1),2)) ', RW: r = ' num2str(round(r(2),2)) ', p = ' num2str(round(p(2),2))]})


subplot(1,2,2)
x = SPeff(refses,:,3);
y = PFeff(refses,:,3);
[r(3),p(3)] = corr(x(:),y(:));
scatter(x(:),y(:),'filled','MarkerFaceColor',[.8 0 0])
hold on
x = SPeff(refses,:,4);
y = PFeff(refses,:,4);
[r(4),p(4)] = corr(x(:),y(:));
scatter(x(:),y(:),'filled','MarkerFaceColor',[0 .8 0])
xlim(xrange)
ylim(yrange)
hold off
xlabel('\Delta Running Speed (cm/s)')
% ylabel('\Delta Place Field Density (AU)')
ylabel('\Delta PV vectors in bin (AU)')
title({['D-D, ' experiment{numexp(refses)}] ['AP: r = ' num2str(round(r(3),2)) ', p = ' num2str(round(p(3),2)) ', RW: r = ' num2str(round(r(4),2)) ', p = ' num2str(round(p(4),2))]})
temppos = get(gca,"position");
legend({'AP', 'RW'},"Location","northeastoutside")
set(gca,"position",temppos)


shift = 0;
readout = 1;% 1 = PF density, 3 = place field centers
xrange = [-1 .5];%[-15 5]
yrange = [-2.5 2];%[-.2 .5];%

PFeff = zeros(6,4,length(numexp));
for i = 1:length(numexp) % Session
    for j = 1:4 % Mouse
        PFeff(1,j,i) = mean(PFVENT(i,int+shift,j,readout))-mean(PFVENT(i,refint,j,readout));
        PFeff(2,j,i) = mean(PFVENT(i,int+10+shift,j,readout))-mean(PFVENT(i,refint,j,readout));
        PFeff(3,j,i) = mean(PFVENT(i,int+20+shift,j,readout))-mean(PFVENT(i,refint,j,readout));
        
        PFeff(4,j,i) = mean(PFDORS(i,int+shift,j,readout))-mean(PFDORS(i,refint,j,readout));
        PFeff(5,j,i) = mean(PFDORS(i,int+10+shift,j,readout))-mean(PFDORS(i,refint,j,readout));
        PFeff(6,j,i) = mean(PFDORS(i,int+20+shift,j,readout))-mean(PFDORS(i,refint,j,readout));
        
    end
end


figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[3 3 [25 15]],...
    'PaperUnits','centimeters',...
    'PaperSize', [25 15])

strcomp = {'I-AP' 'II' 'III-RW'};
yrange = [-.5 .5];
p = zeros(length(numexp),2);
for i = 1:length(numexp)
    refses = i;
    
    subplot(2,length(numexp),i)
    x = PFeff(1:3,:,refses)';
%     [p,tbl,stats] = anova1(x,[],"off");
    [~,p(i,1)] = ttest(x(:,1),x(:,3));
    plot(1:3,x','color',[.8 .8 .8])
    hold on
    boxplot(x,strcomp)
    box('off')
    ylim(yrange)
    title(['I-D, ' experiment{numexp(i)} ', p = ' num2str(round(p(i,1),2))])
    subplot(2,length(numexp),length(numexp)+i)
    x = PFeff(4:6,:,refses)';
%     [p,tbl,stats] = anova1(x,[],"off");
    [~,p(i,2)] = ttest(x(:,1),x(:,3));
    plot(1:3,x','color',[.8 .8 .8])
    hold on
    boxplot(x,strcomp)
    box('off')
    ylim(yrange)
    title(['D-D, ' experiment{numexp(i)} ', p = ' num2str(round(p(i,2),2))])
end

