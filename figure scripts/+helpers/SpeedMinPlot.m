function SpeedMinPlot(CAIM,session)
    
experiment = {'B1' 'B2' 'B3' 'T1' 'T2' 'TN-1' 'TN' 'P1' };
mouseID = {'M259','M261', 'M270','M272','M262','M263','M271','M278'};


%%
numbin = 3; 

numice = 1:8;
counts = zeros(length(numice),length(session),numbin); % 8 mice, 7 sessions
Fracts = zeros(length(numice),length(session),numbin); 

for j = 1:length(session)
  
    for i = 1:length(numice)
     
        speed = CAIM(session(j),numice(i)).behave.speedbin;
        speed = speed(:,1:30);
        
        mn = mean(speed,2, 'omitnan');
        st = nanstd(speed, [], 2);%, 'omitnan');
        z = (speed - mn)./st;

        roundNum = size(z,1);
        roundLen = size(z,2);
        binsize = roundLen/numbin; %mm
        intervals = 1 : binsize :roundLen+1;

        figure(1)
        hold on

        for rnd = 1:roundNum
            new = z(rnd,:);
            new(new>-1)=nan;

            for inv= 1:length(intervals)-1
                [minInterval, minIdx] = nanmin(new(intervals(inv):intervals(inv+1)-1));
                new(intervals(inv):intervals(inv+1)-1) = nan;
                new(intervals(inv)+minIdx-1) = minInterval;  
            end

    %         new(new==100) = nan;
            for k=1:numbin
               counts(i,j,k)= counts(i,j,k) + nansum((new((k-1)*binsize+1 : k*binsize)));
            end 

            plot((1:roundLen)*5, new, 'k.','markersize',10, 'LineWidth',1.5)

        end

        box off
        ylabel('z scored speed')
        xlabel('space bins')
        ylim([-3 0])
        title(['Min. speed zscored'])
        % calculate the n of minimums in bins and normalize 
        Fracts (i,j,:) = (squeeze(counts(i,j,:))/sum(counts(i,j,:)));
    
    end    
end


%% plot selected sessions for selected bin sizes (10 cm , 50 cm...)
figure

Fracts_n = Fracts(:,:,:); 
bar(squeeze(mean(Fracts_n(:,:,1),1)),'facecolor', [0.7 0.7 0.7]) % make a mean across columns (dim=1, mice) and plot
SEM = std(Fracts_n(:,:,1), [],1)/sqrt(size(Fracts_n(:,:,1),1));
hold on
errorbar(squeeze(mean(Fracts_n(:,:,1),1)),SEM)
box off
hold on 
ax = gca;
    ax.XTick = 1:length(session);
    ax.XTickLabel = {experiment{session}};
        ylim([0 1])

ylabel('fraction of min. in the first 50 cm')
title('Fraction of stops 50 cm-airpuff location')
end