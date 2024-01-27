%%
load('/media/2Photon/Negar Nikbakht/ANALYSIS/BigFatCluster.mat')
load('/media/2Photon/Negar Nikbakht/ANALYSIS/BigFatPCA.mat')

pdfpathname = '/media/2Photon/Negar Nikbakht/Final_Plots/';

lnwd = 1;
ftsz = 8;
% experiment = {'Pretrain-1' 'Train-1' 'posttrain-1' 'extinction'};

%%
num = [3]; % sessions

% mean or std?
% For the error the real values are substracted from the fit/prediction.
% This creates a distribution around zero. The std of this distribution is
% used as a measure, how far the prediction is away from the real value.
% For the classifier, since it is binary, the mean of the absolute value of
% this distribition equals the percentage of wrong predictions. 
MoS = 2; 

% wihich interval tolook at: 1:3 learning all/run/rest 4:6 testing all/run/rest 
% For speed and position you would take interval 5
% For the classifier interval 4 (This is included in the loop)
int = 5; 

% Which measures to read out. 1 - sin; 2 - cos; 3 - speed; 4 - classifier; 5 - position
% 6-10 same variables but fitted from PCA traces
Readout = 1:5;

perfVent = [];
perferrVent = [];
perfDors = [];
perferrDors = [];

pSpeed = zeros(length(num),size(CAIM,2));
pClass = zeros(length(num),size(CAIM,2));
pLoc = zeros(length(num),size(CAIM,2));

for i = 1:length(num)
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).A) || (CAIM(num(i),fullCAIM(j)).behave.numrounds)<5
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end

    fullCAIM(emptyCAIM) = [];    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);

        y = PCA(num(i),k).GLMout.mdlperf;
        y = y(Readout,int,MoS);
            y(4,1,1) = PCA(num(i),k).GLMout.mdlperf(4,4,1); % Special line for classifier
        yy = PCA(num(i),k).GLMout.glmshuffleMean;
        yy = yy(Readout,int,MoS);
            yy(4,1,1) = PCA(num(i),k).GLMout.glmshuffleMean(4,4,1); % Special line for classifier
        yyerr = PCA(num(i),k).GLMout.glmshuffleStd;
        yyerr = yyerr(Readout,int,MoS);
            yyerr(4,1,1) = PCA(num(i),k).GLMout.glmshuffleStd(4,4,1); % Special line for classifier
        if k<5              
            perfVent = [perfVent y];
            perferrVent = [perferrVent yy];
        else
            perfDors = [perfDors y];
            perferrDors = [perferrDors yy];
        end
        %%
        glmp = PCA(num(i),k).GLMout.glmp;
        glmp = glmp(Readout,int,MoS);
        glmp(4,1,1) = PCA(num(i),k).GLMout.glmp(4,4,1);
        
        pSpeed(i,j) = glmp(3);
        pClass(i,j) = glmp(4);
        pLoc(i,j) = glmp(5);
    end
end


%%
% figure('color',[1 1 1],...
%     'renderer','painters',...
%     'visible','on',...
%     'Units','centimeters',...
%     'position',[10 2 [ 2*8.9 sqrt(2)*8.9]],...
%     'PaperUnits','centimeters',...
%     'PaperSize', [2*8.9 sqrt(2)*8.9]);
figure 
subplot(1,3,1)

readout = 4;
x = [1 2 4 5];
y = [mean(perfVent(readout,:),2) mean(perferrVent(readout,:),2) mean(perfDors(readout,:),2) mean(perferrDors(readout,:),2)];
yerr = [std(perfVent(readout,:),0,2) std(perferrVent(readout,:),0,2) std(perfDors(readout,:),0,2) std(perferrDors(readout,:),0,2)]/sqrt(4);
b = bar(x,y);
% b.FaceColor = 'flat';
% b.LineStyle = 'none';
% b.CData(1,:) = [.2 .2 .2];
% b.CData(2,:) = [.5 .5 .5];
% b.CData(3,:) = [.2 .2 .2];
% b.CData(4,:) = [.5 .5 .5];

hold on
errorbar(x,y,nan(size(yerr)),yerr,'.',...
            'Marker','none',...
            'Color',[0 0 0],...
            'LineWidth',lnwd)
hold off
box off

ax = gca;
ax.LineWidth = lnwd;
ax.FontSize = ftsz-2;
ax.XLim = [0 3];
% ax.YLim = [0 .015];
ax.XTick = [1.5 4.5];
ax.XTickLabel = {'Ventral' 'Dorsal'};
xlabel('running or resting')
ylabel('% wrong classification')

subplot(1,3,2)
readout = 3;
x = [1 2 4 5];
y = [mean(perfVent(readout,:),2) mean(perferrVent(readout,:),2) mean(perfDors(readout,:),2) mean(perferrDors(readout,:),2)];
yerr = [std(perfVent(readout,:),0,2) std(perferrVent(readout,:),0,2) std(perfDors(readout,:),0,2) std(perferrDors(readout,:),0,2)]/sqrt(4);
b = bar(x,y);
% b.FaceColor = 'flat';
% b.LineStyle = 'none';
% b.CData(1,:) = [.2 .2 .2];
% b.CData(2,:) = [.5 .5 .5];
% b.CData(3,:) = [.2 .2 .2];
% b.CData(4,:) = [.5 .5 .5];

hold on
errorbar(x,y,nan(size(yerr)),yerr,'.',...
            'Marker','none',...
            'Color',[0 0 0],...
            'LineWidth',lnwd)
hold off
box off

ax = gca;
ax.LineWidth = lnwd;
ax.FontSize = ftsz-2;
% ax.XLim = [0 3];
ax.YLim = [0 8];
ax.XTick = [1.5 4.5];
ax.XTickLabel = {'Ventral' 'Dorsal'};

xlabel('speed prediction')
ylabel('std(error) /cm s^{-1}')
title('Baseline Sessions', 'fontsize',ftsz-2)
[Speedh, Speedp] = ttest2 (perfVent(readout,:),perfDors(readout,:)) 

subplot(1,3,3)

%% multiple measure anova
y = [perfVent(readout,:) perferrVent(readout,:)  perfDors(readout,:) perferrDors(readout,:)]';
a = size(perfVent,2);
g1 = {};
g2 = {};
g3 = {};
g1(1:a,1) = {'vent' };
g1(end+1:end+a) =  {'vent' };
g2(1:a,1) = {'data' };
g2(end+1:end+a) =  {'shuf' };
g3(1:4,1)={'day1' };
g3(5:8,1)={'day2' };
g3(9:12,1)={'day3' };
g3 = [g3;g3;g3;g3];
a = size(perfDors,2);
g1(end+1:end+a) =  {'dors' };
g1(end+1:end+a) =  {'dors' };
g2(end+1:end+a) = {'data' };
g2(end+1:end+a) =  {'shuf' };
g3(1:4,1)={'day1' };
g3(5:8,1)={'day2' };
g3(9:12,1)={'day3' };
 
[ASPp,tbl,stats] = anovan (y,{g1 g2,g3},'varnames',{'V/D','Shuf','days'});
% [ASPp,tbl,stats] = anovan (y,{g1 g2,g3},'varnames',{'V/D','Shuf','days'});

[c,m,h,nms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');



%% n way anova 
y = [perfVent(readout,:) perfDors(readout,:)]';
a = size(perfVent,2);
g1 = {};
g2 = {};
g3 = {};
g1(1:a,1) = {'vent' };
g3(1:4,1)={'day1' };
g3(5:8,1)={'day2' };
g3(9:12,1)={'day3' };
g3 = [g3;g3;g3;g3];
a = size(perfDors,2);
g1(end+1:end+a) =  {'dors' };
g3(1:4,1)={'day1' };
g3(5:8,1)={'day2' };
g3(9:12,1)={'day3' };
 
[ASPp,tbl,stats] = anovan (y,{g1 ,g3},'varnames',{'V/D','days'});
% [ASPp,tbl,stats] = anovan (y,{g1 g2,g3},'varnames',{'V/D','Shuf','days'});

[c,m,h,nms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
%% n way anova 
y = [perfVent(readout,:) perfDors(readout,:)]';
a = size(perfVent,2);
g1 = {};
g2 = {};
g3 = {};
g1(1:a,1) = {'vent' };
g3(1:4,1)={'day1' };
g3(5:8,1)={'day2' };
g3(9:12,1)={'day3' };
g3 = [g3;g3;g3;g3];
a = size(perfDors,2);
g1(end+1:end+a) =  {'dors' };
g3(1:4,1)={'day1' };
g3(5:8,1)={'day2' };
g3(9:12,1)={'day3' };
 
[ASPp,tbl,stats] = anovan (y,{g1 ,g3},'varnames',{'V/D','days'});
% [ASPp,tbl,stats] = anovan (y,{g1 g2,g3},'varnames',{'V/D','Shuf','days'});

[c,m,h,nms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');

%%
readout = 5;
y = [mean(perfVent(readout,:),2) mean(perferrVent(readout,:),2) mean(perfDors(readout,:),2) mean(perferrDors(readout,:),2)];
yerr = [std(perfVent(readout,:),0,2) std(perferrVent(readout,:),0,2) std(perfDors(readout,:),0,2) std(perferrDors(readout,:),0,2)]/sqrt(4);
y = 150*y/(pi);
yerr = 150*yerr/(pi);

b = bar(x,y);
% b.FaceColor = 'flat';
% b.LineStyle = 'none';
% b.CData(1,:) = [.2 .2 .2];
% b.CData(2,:) = [.5 .5 .5];
% b.CData(3,:) = [.2 .2 .2];
% b.CData(4,:) = [.5 .5 .5];

hold on
errorbar(x,y,nan(size(yerr)),yerr,'.',...
            'Marker','none',...
            'Color',[0 0 0],...
            'LineWidth',lnwd)
hold off
box off

ax = gca;
ax.LineWidth = lnwd;
ax.FontSize = ftsz-2;
% ax.XLim = [0 3];

% ax.YLim = [0 .015];
ax.XTick = [1.5 4.5];
ax.XTickLabel = {'Ventral' 'Dorsal'};
ylim ([0 14])
xlabel('position prediction')
ylabel('std(error) /cm')

[Posh, Posp] = ttest2 (perfVent(readout,:),perfDors(readout,:)) % ventral pos. vs. dorsal pos. 

[PoshV, PospV] = ttest2 (perfVent(readout,:),(perferrVent(readout,:))) % ventral pos. against shuffle 
[PoshD, PospD] = ttest2 (perfDors(readout,:),(perferrDors(readout,:))) % dorsal pos. against shuffle 


% [ APosp] = anova2 ([perfVent(readout,:); perfDors(readout,:)]',3)
% 
% print(gcf, '-dpdf', [pdfpathname 'GLM_Shuffle_T1.pdf'])
%% multiple measure anova for position
y = [perfVent(readout,:) perferrVent(readout,:)  perfDors(readout,:) perferrDors(readout,:)]';
a = size(perfVent,2);
g1 = {};
g2 = {};
% g3 = {};
g1(1:a,1) = {'vent' };
g1(end+1:end+a) =  {'vent' };
g2(1:a,1) = {'data' };
g2(end+1:end+a) =  {'shuf' };
% g3(1:4,1)={'day1' };
% g3(5:8,1)={'day2' };
% g3(9:12,1)={'day3' };
g3 = [g3;g3;g3;g3];
a = size(perfDors,2);
g1(end+1:end+a) =  {'dors' };
g1(end+1:end+a) =  {'dors' };
g2(end+1:end+a) = {'data' };
g2(end+1:end+a) =  {'shuf' };
 
[APosp,tbl,stats] = anovan (y,{g1 g2},'varnames',{'V/D','Shuf'});
[c,m,h,nms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
%% n way anova 
y = [perfVent(readout,:) perfDors(readout,:)]';
a = size(perfVent,2);
g1 = {};
g2 = {};
g3 = {};
g1(1:a,1) = {'vent' };
g3(1:4,1)={'day1' };
g3(5:8,1)={'day2' };
g3(9:12,1)={'day3' };
g3(9:12,1)={'day4' };
g3 = [g3;g3;g3;g3];
a = size(perfDors,2);
g1(end+1:end+a) =  {'dors' };
g3(1:4,1)={'day1' };
g3(5:8,1)={'day2' };
g3(9:12,1)={'day3' };
 
[ASPp,tbl,stats] = anovan (y,{g1 ,g3},'varnames',{'V/D','days'});
% [ASPp,tbl,stats] = anovan (y,{g1 g2,g3},'varnames',{'V/D','Shuf','days'});

[c,m,h,nms] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');

%% p-vlaues of linear model learning

num = [1:3]; % sessions
MoS = 2; % mean or std
int = 5; % wihich interval to to test 1:3 learning all/run/rest 4:6 testing all/run/rest
% Which measures to read out. 1 - sin; 2 - cos; 3 - speed; 4 - classifier; 5 - position
% 6-10 same variables but fitted from PCA traces
Readout = 1:5;

glmp = NaN(5,length(num),2,size(CAIM,2));
mdlperf = NaN(5,length(num),2,size(CAIM,2));
numrounds = NaN(length(num),size(CAIM,2));
numcells = NaN(length(num),size(CAIM,2));
compperf = zeros(5,20);


for i = 1:length(num)
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(PCA(num(i),fullCAIM(j)).GLMout)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        numrounds(i,k) = CAIM(num(i),k).behave.numrounds; 
        numcells(i,k) = size(CAIM(num(i),k).A,2);
        GLMout = PCA(num(i),k).GLMout;
        if ~isempty(GLMout)
            
            y = GLMout.glmp;
            % 1.var 2.ExpID 3.mean or std 4.mouseID
            glmp(:,i,:,k) = y(Readout,int,:); 
            glmp(4,i,:,k) = y(Readout(4),int-1,:);                      
        end
    end
end

%% p-value plot
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

subplot(3,3,1)
y = permute(glmp(1,:,2,1:4),[4 2 1 3]);
boxplot(y)
title(['p of sin, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
ax = gca;
% ax.XTickLabel = experiment(3:6);

subplot(3,3,2)
y = permute(glmp(2,:,2,1:4),[4 2 1 3]);
boxplot(y)
title(['p of cos, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
ax = gca;
% ax.XTickLabel = experiment(3:6);

subplot(3,3,3)
y = permute(glmp(5,:,2,1:4),[4 2 1 3]);
boxplot(y)
title(['p of theta, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
y
ax = gca;
% ax.XTickLabel = experiment(3:6);

subplot(3,2,3)
y = permute(glmp(3,:,2,1:4),[4 2 1 3]);
boxplot(y)
title(['p of speed, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
y
ax = gca;
% ax.XTickLabel = experiment(3:6);

subplot(3,2,4)
y = permute(glmp(4,:,1,1:4),[4 2 1 3]);
boxplot(y)
title(['p of classifier, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
ax = gca;
% ax.XTickLabel = experiment(3:6);
% 
% subplot(3,2,5)
% y = permute(glmp(5,:,2,1:4),[2 4 1 3]);
% scatter(numrounds(:),y(:))
% ax = gca;
% xlabel('rounds/session')

% subplot(3,2,6)
% y = permute(glmp(5,:,2,1:4),[2 4 1 3]);
% scatter(numcells(:),y(:))
% ax = gca;
% xlabel('cells/session')

print(gcf, '-dpdf','GLM p-value2')
% print(gcf, '-dpdf', [pdfpathname 'GLM p-value2.pdf'])
%% P-values for dorsal
filename = '/media/2Photon/Negar Nikbakht/Final_Plots/Stats/GLMpValues.xlsx';
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

subplot(3,3,1)
y = permute(glmp(1,:,2,5:8),[4 2 1 3]);
boxplot(y)
title(['p of sin, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
ax = gca;
% ax.XTickLabel = experiment(3:6);

subplot(3,3,2)
y = permute(glmp(2,:,2,5:8),[4 2 1 3]);
boxplot(y)
title(['p of cos, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
ax = gca;
% ax.XTickLabel = experiment(3:6);

subplot(3,3,3)
y = permute(glmp(5,:,2,5:8),[4 2 1 3]);
boxplot(y)
title(['p of theta, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% writematrix(y,filename,'Sheet',1,'Range', 'b2:d5')
% ylim([-.05 1.05]);
y
ax = gca;
% ax.XTickLabel = experiment(3:6);

subplot(3,2,3)
y = permute(glmp(3,:,2,5:8),[4 2 1 3]);
boxplot(y)
title(['p of speed, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
y
ax = gca;
% ax.XTickLabel = experiment(3:6);

subplot(3,2,4)
y = permute(glmp(4,:,1,5:8),[4 2 1 3]);
boxplot(y)
title(['p of classifier, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
% ylim([-.05 1.05]);
ax = gca;
% ax.XTickLabel = experiment(3:6);
% 
% subplot(3,2,5)
% y = permute(glmp(5,:,2,1:4),[2 4 1 3]);
% scatter(numrounds(:),y(:))
% ax = gca;
% xlabel('rounds/session')

% subplot(3,2,6)
% y = permute(glmp(5,:,2,5:8),[2 4 1 3]);
% scatter(numcells(:),y(:))
% ax = gca;
% xlabel('cells/session')

print(gcf, '-dpdf','GLM p-value2')
% print(gcf, '-dpdf', [pdfpathname 'GLM p-value2.pdf'])
%% individual component number and shuffle

int = 5;
% 1.var 2.ExpID 3.mean/std 4.mouseID
Readout = 1:5;

num = [1];
for i = 1:length(num)
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).network) || isempty(CAIM(num(i),fullCAIM(j)).plcfield)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end

    fullCAIM(emptyCAIM) = [];    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);

        x = PCA(num(i),k).GLMout.numcomp;
        y = PCA(num(i),k).GLMout.compperf;
        MoS = 2;
        y = y(Readout,int,MoS,:);
             y(4,1,1,:) = PCA(num(i),k).GLMout.compperf(4,4,1,:);
        y = permute(y,[1 4 2 3]);
        yy = PCA(num(i),k).GLMout.compshuffleMean;
        yy = yy(Readout,int,MoS,:,:);
            yy(4,1,1,:,:) = PCA(num(i),k).GLMout.compshuffleMean(4,4,1,:);
        yyerr = std(yy,0,5);
        yyerr = permute(yyerr,[1 4 2 3]);
%         yy = mean(yy,5);
        yy = permute(yy,[1 4 2 3]);

        figure('color',[1 1 1],...
                'renderer','painters',...
                'visible','on',...
                'Units','centimeters',...
                'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
                'PaperUnits','centimeters',...
                'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
        for ii = 1:length(Readout)            
            subplot(5,1,ii)
            hold on
%             axes('position',[9 9 5  5])

            plot(x,y(ii,:))
            errorbar(x,yy(ii,:),yyerr(ii,:))
        end
        subplot(5,1,1)
        title(['sin position'])
        ylabel('std(predict(sin)-real(sin))')
        xlabel('number of components')
        subplot(5,1,2)
        title('cos position')
        ylabel('std(predict(cos)-real(cos))')
        xlabel('number of components')
        subplot(5,1,3)
        title('Speed')
        ylabel('std(predict(speed)-real(speed))')
        xlabel('number of components')
        subplot(5,1,4)
        title('Predictor')
        ylabel('% false predictions')
        xlabel('number of components')
        subplot(5,1,5)
        title('position')
%         axes('position',[900 900 5  5])
        box on
        ylabel('std(predict(ang)-real(ang))')
        xlabel('number of components')
        print(gcf, '-dpdf',['GLM performance 2' mouse{k}])
    end
end

