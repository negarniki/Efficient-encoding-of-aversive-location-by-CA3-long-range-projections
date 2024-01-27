function PlaceCoderPlot(CAIM,cclustID,session,numice)
% In this code we calculate the fraction of place coding axons to all fibers in
% % the given session and make a bar plot

experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'Tn-1' 'Tn' 'P1'};

numPlcD = zeros(length(session),size(CAIM,2)-4,2);
numPlcV = zeros(length(session),size(CAIM,2)-4,2);
% This loop goes through the rows (the experiemts)
for i = 1:length(session)   
    for j = 1:length(numice)    
        k = numice(j);  
        cclust = CAIM(session(i),k).cclust;
        % losonczy criterium
%         isplace = cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcvctp)<=.05;
        % Dombeck criterium
        isplace = cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcfldp)<=.05;
        % all cells 
%         isplace = cclust(:,cclustID.plcvct)>0%& cclust(:,cclustID.plcvctp)>=.05;
        
        if j < 5
            numPlcV(i,j,:) = [sum(isplace) length(isplace)];
        else
            numPlcD(i,j-4,:) = [sum(isplace) length(isplace)];
        end
    end
end

%%
figSize=[5 5];
figure('color',[1 1 1],...
      'renderer','painters',...
      'visible','on',...
      'Units','centimeters',...
      'position',[20 5 figSize ],...
      'PaperUnits','centimeters',...
      'PaperSize', figSize )

% plot fraction/percentage of place coding axons to the entire number of components/axons

a = sum(numPlcV(:,:,1),2); 
b = sum(numPlcD(:,:,1),2);
c = sum(numPlcV(:,:,2),2);
d = sum(numPlcD (:,:,2),2);

fPlcV = a./c;
fPlcD = b./d;

yV = [];
yD = [];
for i = 1:length(session)
    yV = [yV mean(fPlcV(i,:))*100];
    yD = [yD mean(fPlcD(i,:))*100];
end
bar([1:length(session)], yV,0.5,'FaceColor',[0 1 0]) % choose which sessions to plot
hold on
bar([1:length(session)]+length(session), yD,0.5,'FaceColor',[1 1 0]) 

box off 
ax = gca;
ax.XTick = 1:9;
% ax.XTick = 1:3;
xtickangle(23)
ax.XTickLabel = {experiment{session}}; 

% ylim([0 1])
ylabel ('% of significant place coders')
title ('Fraction of Place Coders ')


