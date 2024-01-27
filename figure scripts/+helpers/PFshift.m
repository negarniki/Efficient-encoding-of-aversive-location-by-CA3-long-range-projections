function [plcVentral, plcDorsal, negarD, negarV,DNrestActFract, pooledf] = PFshift(CAIMcorr,cclustID,numses)
% this code calculates theplace cell stabillity and how cells have shifted
% from one session to another 

% pdfpathname = '/media/2Photon/Negar Nikbakht/Final_Plots/PFShift/ConsecutiveSess/';
% load('/media/2Photon/Negar Nikbakht/ANALYSIS/cclustID','cclustID')

%% Shift of center of mass of PCs
% sessions. First session for PC definition. second session for the shift
% numses = [1 2];
% all mice are considered
numice = [1:4; 5:8];
for j = 1:2
    % Since we want a subset of mice we make a new cclust
    cclust = [];
    cclustref = [];
    for i = numice(j,:)
        
        tempref = CAIMcorr(numses(1),i).cclust;
        temp =  CAIMcorr(numses(2),i).cclust;
        temp = temp(1:size(tempref,1),:);
        cclustref = [cclustref; tempref];
        cclust = [cclust; temp];
        
    end

    % define pcs
    isplace = cclustref(:,cclustID.plcfld)>0 & cclustref(:,cclustID.plcfldp)<.05;
    isplace2 = cclust(:,cclustID.plcfld)>0 & cclust(:,cclustID.plcfldp)<.05;
    isthere =  cclust(:,cclustID.expID)~=0 & ~isnan(cclust(:,cclustID.expID)) & cclustref(:,cclustID.expID)~=0 & ~isnan(cclustref(:,cclustID.expID));
    isDeNovo = isthere & ~isplace & isplace2;
    isStable = isthere & isplace & isplace2;
    isDisap = isthere & ~isDeNovo & ~isStable;
    
    plccenter = cclustref(:,cclustID.plcfld);
    plccenter(:,2) = cclust(:,cclustID.plcfld);       
    
    % This was missing before, so all cells were included
    plccenter = plccenter(isplace&isplace2,:);
    
    plccenter = plccenter(:,1) - plccenter(:,2);
    plccenter(isnan(plccenter)) = [];
    plccenter(plccenter<-750) = plccenter(plccenter<-750)+1500;   
    plccenter(plccenter>750) = plccenter(plccenter>750)-1500;
    plccenter = round((plccenter)/10);
%     plccenter = abs(plccenter);
    
    % read out frequencies
    meanf = cclust(isthere,cclustID.meanf+(0:2));
    meanfref = cclustref(isthere,cclustID.meanf+(0:2));
    meanf(meanf==0 & meanfref == 0) = nan;
    meanfPC = cclust(isStable,cclustID.meanf+(0:2));
    meanfPCref = cclustref(isStable,cclustID.meanf+(0:2));
    meanfPC(meanfPC==0 & meanfPCref == 0) = nan;
    meanfDN = cclust(isDeNovo,cclustID.meanf+(0:2));
    meanfDNref = cclustref(isDeNovo,cclustID.meanf+(0:2));
    meanfDN(meanfDN==0 & meanfDNref == 0) = nan;
    meanfDS = cclust(isDisap,cclustID.meanf+(0:2));
    meanfDSref = cclustref(isDisap,cclustID.meanf+(0:2));
    meanfDS(meanfDS==0 & meanfDSref == 0) = nan;

    % We check for a pure rest active, 10 % pref or general pref
    meanfrunall = (cclustref(:,cclustID.meanf+1)./cclustref(:,cclustID.meanf+2)); %Calculate quotient
    meanfrunall = meanfrunall(isDeNovo);
    
    if j == 1
        plcVentral = plccenter;
        fDivVentral = meanf - meanfref;
        fPCDivVentral = meanfPC - meanfPCref;
        fDNDivVentral = meanfDN - meanfDNref;
        fDSDivVentral = meanfDS - meanfDSref;
        fVentral = meanf;
        fPCVentral = meanfPC;
        fDNVentral = meanfDN;
        fDSVentral = meanfDS;
        frefVentral =  meanfref;
        fPCrefVentral =  meanfPCref;
        fDNrefVentral =  meanfDNref;
        fDSrefVentral =  meanfDSref;
        % Here the logical is made for denovo and rest prefenernce
        DNrestActVentral = [meanfrunall == 0 meanfrunall<.1 meanfrunall < 1]; 
    else
        plcDorsal = plccenter;
        fDivDorsal = meanf - meanfref;
        fPCDivDorsal = meanfPC - meanfPCref; 
        fDNDivDorsal = meanfDN - meanfDNref; 
        fDSDivDorsal = meanfDS - meanfDSref; 
        fDorsal = meanf;
        fPCDorsal = meanfPC;
        fDNDorsal = meanfDN;
        fDSDorsal = meanfDS;
        frefDorsal =  meanfref;
        fPCrefDorsal =  meanfPCref;
        fDNrefDorsal =  meanfDNref;
        fDSrefDorsal =  meanfDSref;
        DNrestActDorsal = [meanfrunall == 0 meanfrunall<.1 meanfrunall < 1];
    end
    
    % read out place preference vectors   
end

DNrestActFract = [nansum(DNrestActVentral)' nansum(DNrestActDorsal)';
    size(DNrestActVentral,1) size(DNrestActDorsal,1)]

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[15 1 [18.3 26]],...
    'PaperUnits','centimeters',...
    'PaperSize', [18.3 26])

binSize = 5;
bins = -75-binSize/2:binSize:75+binSize/2;

subplot(5,2,1)
hold off
histogram(plcVentral,bins);%,'normalization','probability');
hold on
histogram(plcDorsal,bins);%,'normalization','probability');
xlim([bins(1) bins(end)])
xlabel('PF shift (cm)')
ylabel('counts')
legend({'Ventral' 'Dorsal'},'location','northeast')
title('PF shift')

subplot(5,2,3)
plcVcount = histcounts(plcVentral,bins,'normalization','probability');
plcDcount = histcounts(plcDorsal,bins,'normalization','probability');
hold off
x = bins(1:end-1)+binSize/2;
plot(x,cumsum(plcVcount))
hold on
plot(x,cumsum(plcDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('PF shift (cm)')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
title('PF shift probability')
grid on

% Ststistical numbers:
StatSign = [mean(plcVentral) std(plcVentral) median(plcVentral);
mean(plcDorsal) std(plcDorsal) median(plcDorsal)];


bins = -binSize/2:binSize:75+binSize/2;
x = bins(1:end-1)+binSize/2;

subplot(5,2,2)
hold off
histogram(abs(plcVentral),bins);%,'normalization','probability');
hold on
histogram(abs(plcDorsal),bins);%,'normalization','probability');
xlim([bins(1) x(end)])
xlabel('PF shift (cm)')
ylabel('counts')
legend({'Ventral' 'Dorsal'},'location','northeast')
title('PF shift')

subplot(5,2,4)
plcVcount = histcounts(abs(plcVentral),bins,'normalization','probability');
plcVcount = cumsum(plcVcount);
plcDcount = histcounts(abs(plcDorsal),bins,'normalization','probability');
plcDcount = cumsum(plcDcount);
hold off

plot(x,(plcVcount))
hold on
plot(x,(plcDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('PF shift (cm)')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
title('PF shift probability')
grid on

%
bins = -5:.02:5;

subplot(5,3,7)
x = fDivVentral(:,1);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fDivDorsal(:,1);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('f-Change')
grid on

subplot(5,3,8)
x = fDivVentral(:,2);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fDivDorsal(:,2);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('f-Change - run')
grid on

subplot(5,3,9)
x = fDivVentral(:,3);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fDivDorsal(:,3);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('f-Change - rest')
grid on

subplot(5,3,10)
x = fPCDivVentral(:,1);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCDivDorsal(:,1);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('PC f-Change')
grid on

subplot(5,3,11)
x = fPCDivVentral(:,2);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCDivDorsal(:,2);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('PC f-Change - run')
grid on

subplot(5,3,12)
x = fPCDivVentral(:,3);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCDivDorsal(:,3);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('PC f-Change - rest')
grid on

subplot(5,3,13)
x = fDNDivVentral(:,1);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fDNDivDorsal(:,1);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('DeNovo f-Change')
grid on

subplot(5,3,14)
x = fDNDivVentral(:,2);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fDNDivDorsal(:,2);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('DeNovo f-Change - run')
grid on

subplot(5,3,15)
x = fDNDivVentral(:,3);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fDNDivDorsal(:,3);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount))
hold on
plot(x,(fDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'Ventral' 'Dorsal'},'location','southeast')
legend('boxoff')
title('DN f-Change - rest')
grid on

%% comparissons within groups
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[15 1 [18.3 26]],...
    'PaperUnits','centimeters',...
    'PaperSize', [18.3 26])

bins = -6:.02:6;
subplot(6,3,1)
x = fDivVentral(:,1);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCDivVentral(:,1);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNDivVentral(:,1);
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSDivVentral(:,1);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'All' 'PCs' 'DNs' 'DSs'},'location','southeast')
legend('boxoff')
title('Ventral f-Change')
grid on

subplot(6,3,2)
x = fDivVentral(:,2);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCDivVentral(:,2);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNDivVentral(:,2);
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSDivVentral(:,2);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Ventral running-f-Change')
grid on

subplot(6,3,3)

x = fDivVentral(:,3);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCDivVentral(:,3);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNDivVentral(:,3);
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSDivVentral(:,3);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Ventral rest f-Change')
grid on

subplot(6,3,4)
x = fDivDorsal(:,1);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCDivDorsal(:,1);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNDivDorsal(:,1);
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSDivDorsal(:,1);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
legend('boxoff')
title('Dorsal f-Change')
grid on

subplot(6,3,5)
x = fDivDorsal(:,2);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCDivDorsal(:,2);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNDivDorsal(:,2);
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSDivDorsal(:,2);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Dorsal run f-Change')
grid on

subplot(6,3,6)
x = fDivDorsal(:,3);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCDivDorsal(:,3);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNDivDorsal(:,3);
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSDivDorsal(:,3);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('\Delta f')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Dorsal rest-f-Change')
grid on

% frequencies day 1
bins = 0:.02:10;
subplot(6,3,7)
x = frefVentral(:,1);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCrefVentral(:,1);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNrefVentral(:,1);
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSrefVentral(:,1);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
legend('boxoff')
title('Ventral f, day 1')
grid on
negarV = zeros(1,1);
subplot(6,3,8)
x = frefVentral(:,2);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCrefVentral(:,2);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNrefVentral(:,2);
negarV(1) = sum(x==0); % negar added 
negarV(2)= sum(~isnan(x)); %negar added
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSrefVentral(:,2);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Ventral f, run, day 1')
grid on
%
subplot(6,3,9)
x = frefVentral(:,3);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCrefVentral(:,3);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNrefVentral(:,3);
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSrefVentral(:,3);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
%  fPCVRest=fPCVcount ; % Negar added - to pool the frequencies 
%  fDNVRest= fDNVcount;
% fDSVrest = fDSVcount;

title('Ventral f, rest, day 1')
grid on


subplot(6,3,10)
x = frefDorsal(:,1);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCrefDorsal(:,1);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNrefDorsal(:,1);
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSrefDorsal(:,1);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
legend('boxoff')
title('Dorsal f, day 1')
grid on

negarD = zeros(1,1);
subplot(6,3,11)
x = frefDorsal(:,2);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCrefDorsal(:,2);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNrefDorsal(:,2);
negarD(1) = sum(x==0); % negar added 
negarD(2)= sum(~isnan(x)); %negar added
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSrefDorsal(:,2);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Dorsal f, run, day 1')
grid on


subplot(6,3,12)
x = frefDorsal(:,3);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCrefDorsal(:,3);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNrefDorsal(:,3);
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSrefDorsal(:,3);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Dorsal f, rest, day1')
grid on

% frequencies day 2

% bins = 0:.02:5;
subplot(6,3,13)
x = fVentral(:,1);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCVentral(:,1);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNVentral(:,1);
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSVentral(:,1);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
legend('boxoff')
title('Ventral f, day 2')
grid on

subplot(6,3,14)
x = fVentral(:,2);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCVentral(:,2);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNVentral(:,2);
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSVentral(:,2);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Ventral f, run, day 2')
grid on

subplot(6,3,15)
x = fVentral(:,3);
fVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fVcount = cumsum(fVcount);
x = fPCVentral(:,3);
fPCVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCVcount = cumsum(fPCVcount);
x = fDNVentral(:,3);
fDNVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNVcount = cumsum(fDNVcount);
x = fDSVentral(:,3);
fDSVcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSVcount = cumsum(fDSVcount);

x = bins(1:end-1)+.01;
plot(x,(fVcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCVcount))
plot(x,(fDNVcount))
plot(x,(fDSVcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Ventral f, rest, day 2')
grid on

subplot(6,3,16)
x = fDorsal(:,1);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCDorsal(:,1);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNDorsal(:,1);
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSDorsal(:,1);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
legend('boxoff')
title('Dorsal f, day 2')
grid on

subplot(6,3,17)
x = fDorsal(:,2);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCDorsal(:,2);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNDorsal(:,2);
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSDorsal(:,2);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Dorsal f, run, day 2')
grid on


subplot(6,3,18)
x = fDorsal(:,3);
fDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDcount = cumsum(fDcount);
x = fPCDorsal(:,3);
fPCDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fPCDcount = cumsum(fPCDcount);
x = fDNDorsal(:,3);
fDNDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDNDcount = cumsum(fDNDcount);
x = fDSDorsal(:,3);
fDSDcount = histcounts(x(~isnan(x)),bins,'normalization','probability');
fDSDcount = cumsum(fDSDcount);

x = bins(1:end-1)+.01;
plot(x,(fDcount),'color',[.5 .5 .5])
hold on
plot(x,(fPCDcount))
plot(x,(fDNDcount))
plot(x,(fDSDcount))
xlim([x(1) x(end)])
ylim([0 1])
xlabel('f (events/min)')
ylabel('cum. prob.')
% legend({'All' 'PCs' 'DNs' 'Other'},'location','southeast')
% legend('boxoff')
title('Dorsal f, rest, day 2')
grid on
%% Statistical numbers:
StatAbs = [mean(abs(plcVentral)) std(abs(plcVentral)) median(abs(plcVentral));
mean(abs(plcDorsal)) std(abs(plcDorsal)) median(abs(plcDorsal))];

StatPerc =[x(find(plcVcount>.5,1)) x(find(plcVcount>.66,1)) x(find(plcVcount>.95,1));
x(find(plcDcount>.5,1)) x(find(plcDcount>.66,1)) x(find(plcDcount>.95,1))];

%% pooled output
% all cells Dorsal
pooledf.fDivDorsal = fDivDorsal; % Difference of day2-day1
pooledf.fDorsal= fDorsal; % frequencies day 2
pooledf.frefDorsal= frefDorsal; % frequencies day 1

% PCs Dorsal
pooledf.fPCDivDorsal = fPCDivDorsal; % Difference of day2-day1
pooledf.fPCDorsal= fPCDorsal; % frequencies day 2
pooledf.fPCrefDorsal= fPCrefDorsal; % frequencies day 1

% de novo Dorsal
pooledf.fDNDivDorsal = fDNDivDorsal; % Difference of day2-day1
pooledf.fDNDorsal= fDNDorsal; % frequencies day 2
pooledf.fDNrefDorsal= fDNrefDorsal; % frequencies day 1

% Others Dorsal
pooledf.fDSDivDorsal = fDSDivDorsal; % Difference of day2-day1
pooledf.fDSDorsal= fDSDorsal; % frequencies day 2
pooledf.fDSrefDorsal= fDSrefDorsal; % frequencies day 1

% all Ventral
pooledf.fDivVentral = fDivVentral; % Difference of day2-day1
pooledf.fVentral = fVentral; % frequencies day 2
pooledf.frefVentral = frefVentral; % frequencies day 1

% PCs ventral
pooledf.fPCDivVentral = fPCDivVentral; % Difference of day2-day1
pooledf.fPCVentral = fPCVentral; % frequencies day 2
pooledf.fPCrefVentral = fPCrefVentral; % frequencies day 1

% De novo ventral
pooledf.fDNDivVentral = fDNDivVentral; % Difference of day2-day1
pooledf.fDNVentral = fDNVentral; % frequencies day 2
pooledf.fDNrefVentral = fDNrefVentral; % frequencies day 1

% Others ventral
pooledf.fDSDivVentral = fDSDivVentral; % Difference of day2-day1
pooledf.fDSVentral = fDSVentral; % frequencies day 2
pooledf.fDSrefVentral = fDSrefVentral; % frequencies day 1
end