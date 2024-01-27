% Figure 5

load('BigFatCluster.mat')
load('BigFatPCA.mat')
load('cclustID.mat')

experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'Tn-1' 'Tn' 'P1'};
mouseID = {'M259','M261', 'M270','M272','M262','M263','M271','M278'};

%% A,B
session = [3 4 7 8];
numice = [1 5];
figsize = [20 5];
        
for i = 1:length(session)
    for j = 1:length(numice)
        figure('color',[1 1 1],...
            'renderer','painters',...
            'visible','on',...
            'Units','centimeters',...
            'position',[10 2 figsize],...
            'PaperUnits','centimeters',...
            'PaperSize', figsize);
        ts = CAIM(session(i),numice(j)).behave.tsscn/1000;
        movement = PCA(session(i),numice(j)).GLMout.movement;
        pred = PCA(session(i),numice(j)).GLMout.pred;
        y = 150*(movement(:,5)+pi)/(2*pi);
        plot(ts,y,'color',[0 1 1])
        hold on
        y = 150*(pred(:,10)+pi)/(2*pi);
        plot(ts,y,'color',[0 0 0])
        ylabel('position (cm)')
        xlabel('time (s)')
        title(['mouse ' mouseID{j} ', session ' experiment{session(i)}])
        box off
    end
end

%% C

session = [1:3]; %for 3 baselines 

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

pSpeed = zeros(length(session),size(CAIM,2));
pClass = zeros(length(session),size(CAIM,2));
pLoc = zeros(length(session),size(CAIM,2));

for i = 1:length(session)
     
    for j = 1:size(CAIM,2)    

        y = PCA(session(i),j).GLMout.mdlperf;
        y = y(Readout,int,MoS);
            y(4,1,1) = PCA(session(i),j).GLMout.mdlperf(4,4,1); % Special line for classifier
        yy = PCA(session(i),j).GLMout.glmshuffleMean;
        yy = yy(Readout,int,MoS);
            yy(4,1,1) = PCA(session(i),j).GLMout.glmshuffleMean(4,4,1); % Special line for classifier
        yyerr = PCA(session(i),j).GLMout.glmshuffleStd;
        yyerr = yyerr(Readout,int,MoS);
            yyerr(4,1,1) = PCA(session(i),j).GLMout.glmshuffleStd(4,4,1); % Special line for classifier
        if j<5              
            perfVent = [perfVent y];
            perferrVent = [perferrVent yy];
        else
            perfDors = [perfDors y];
            perferrDors = [perferrDors yy];
        end
        %%
        glmp = PCA(session(i),j).GLMout.glmp;
        glmp = glmp(Readout,int,MoS);
        glmp(4,1,1) = PCA(session(i),j).GLMout.glmp(4,4,1);
        
        pSpeed(i,j) = glmp(3);
        pClass(i,j) = glmp(4);
        pLoc(i,j) = glmp(5);
    end
end

figsize = [7 10];
lnwd = 1;
ftsz = 8;
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 figsize],...
    'PaperUnits','centimeters',...
    'PaperSize', figsize);

readout = 5;
x = [1 2 4 5];
y = [mean(perfVent(readout,:),2) mean(perferrVent(readout,:),2) mean(perfDors(readout,:),2) mean(perferrDors(readout,:),2)];
yerr = [std(perfVent(readout,:),0,2) std(perferrVent(readout,:),0,2) std(perfDors(readout,:),0,2) std(perferrDors(readout,:),0,2)]/sqrt(4);
y = 150*y/(pi);
yerr = 150*yerr/(pi);

b = bar(x,y);
b.FaceColor = 'flat';
b.LineStyle = 'none';
b.CData(1,:) = [.2 .7 .2];
b.CData(2,:) = [.5 .5 .5];
b.CData(3,:) = [.7 .7 .2];
b.CData(4,:) = [.5 .5 .5];

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
ax.XTick = [1.5 4.5];
ax.XTickLabel = {'Ventral' 'Dorsal'};
ylim ([0 12])
xlabel('position prediction')
ylabel('std(error) /cm')

%% D

session = [3 4 7 8]; %for 3 baselines 

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

pSpeed = zeros(length(session),size(CAIM,2));
pClass = zeros(length(session),size(CAIM,2));
pLoc = zeros(length(session),size(CAIM,2));

for i = 1:length(session)
     
    for j = 1:size(CAIM,2)    

        y = PCA(session(i),j).GLMout.mdlperf;
        y = y(Readout,int,MoS);
            y(4,1,1) = PCA(session(i),j).GLMout.mdlperf(4,4,1); % Special line for classifier
        yy = PCA(session(i),j).GLMout.glmshuffleMean;
        yy = yy(Readout,int,MoS);
            yy(4,1,1) = PCA(session(i),j).GLMout.glmshuffleMean(4,4,1); % Special line for classifier
        yyerr = PCA(session(i),j).GLMout.glmshuffleStd;
        yyerr = yyerr(Readout,int,MoS);
            yyerr(4,1,1) = PCA(session(i),j).GLMout.glmshuffleStd(4,4,1); % Special line for classifier
        if j<5              
            perfVent = [perfVent y];
            perferrVent = [perferrVent yy];
        else
            perfDors = [perfDors y];
            perferrDors = [perferrDors yy];
        end
        %%
        glmp = PCA(session(i),j).GLMout.glmp;
        glmp = glmp(Readout,int,MoS);
        glmp(4,1,1) = PCA(session(i),j).GLMout.glmp(4,4,1);
        
        pSpeed(i,j) = glmp(3);
        pClass(i,j) = glmp(4);
        pLoc(i,j) = glmp(5);
    end
end

figsize = [7 7];
lnwd = 1;
ftsz = 8;
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 figsize],...
    'PaperUnits','centimeters',...
    'PaperSize', figsize);

readout = 5;
x = [1:4 6:9];
y = [mean(perfVent(readout,1:4),2)  mean(perfVent(readout,5:8),2) mean(perfVent(readout,9:12),2) mean(perfVent(readout,13:16),2) mean(perfDors(readout,1:4),2) mean(perfDors(readout,5:8),2) mean(perfDors(readout,9:12),2) mean(perfDors(readout,13:16),2)];
yerr =  [std(perfVent(readout,1:4),0,2)  std(perfVent(readout,5:8),0,2) std(perfVent(readout,9:12),0,2) std(perfVent(readout,13:16),0,2) std(perfDors(readout,1:4),0,2) std(perfDors(readout,5:8),0,2) std(perfDors(readout,9:12),0,2) std(perfDors(readout,13:16),0,2) ]/sqrt(4);
y = 150*y/(pi);
yerr = 150*yerr/(pi);

b = bar(x,y);
b.FaceColor = 'flat';
b.LineStyle = 'none';
for i = 1:4
    b.CData(i,:) = [.2 .7 .2];
    b.CData(i+4,:) = [.7 .7 .2];
end
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
ax.XTick = x;
ax.XTickLabel = {experiment{session}};
ylim ([0 12])
xlabel('position prediction')
ylabel('std(error) /cm')



