function [PlacePrefVect,R] = TrackedFeaturePlot(CAIMcorr,cclustID,numexp,DoPlot)
%%
% Output: PlacePrefVect -   columns correspont to the four vectors:
%                           Stable-Ventral DeNovo_Ventral Stable_Dorsal _DeNovo_Dorsal
%                           rows  correspond to read outs
%                           1. length 
%                           2. Std of legnth 
%                           3. SEM of length 
%                           4. angular distance to AP (in degree) 
%                           5. Std of angle 
%                           6. SEM of angle
%                           7. Fractions of feature coders
%                           Third dimension corresponds to transition
%           R - Correlations between PF density and Speed
%                           Dim1 days and mean, 
%                           dim2 I-D vs D-D, 
%                           dim3 re-ocurring vs newly coding
%%
numbin = 1500;
Alpha = linspace(pi,-pi,numbin*1);
DZ = exp(Alpha(numbin/3)*1j);
Z2 = exp(Alpha(numbin/3*2)*1j);
Z3 = exp(Alpha(numbin)*1j);
PF = zeros(length(numexp),size(CAIMcorr,2));

x = 2.5:150/30:150;

% Dim1 days and mean, dim2 I-D vs D-D, dim3 re-ocurring vs newly coding
R = zeros(5,2,2);

% Add up numbers of feature coders to calculate fractions
Nv = zeros(3,size(CAIMcorr,1)-1);
Nd = zeros(3,size(CAIMcorr,1)-1);

for i = 1:length(numexp)-1
    PVvStbl = [];
    PVvRemap = [];
    PVvDeno = [];
    PVdStbl = [];
    PVdRemap = [];
    PVdDeno = [];
    
    PFvStbl = [];
    PFvRemap = [];
    PFvDeno = [];
    PFdStbl = [];
    PFdRemap = [];
    PFdDeno = [];
    
    speedbinV = [];
    speedbinD = [];
    

    for j = 1:size(CAIMcorr,2)
        %% define groups of tracked fibers 
        % i =1;j=1;
        
        cclustref = CAIMcorr(numexp(i),j).cclust;
        cclust = CAIMcorr(numexp(i+1),j).cclust;
        cclust = cclust(1:size(cclustref,1),:);
        isthere =  cclust(:,cclustID.expID)~=0 & ~isnan(cclust(:,cclustID.expID)) & cclustref(:,cclustID.expID)~=0 & ~isnan(cclustref(:,cclustID.expID));
        
        isplace = cclust(:,cclustID.plcfld)>0 & cclust(:,cclustID.plcfldp)<=.05;
        isplaceref = cclustref(:,cclustID.plcfld)>0 & cclustref(:,cclustID.plcfldp)<=.05;
        
        plccenter = cclust(:,cclustID.plcfld);
        plccenter = round(plccenter);
        plccenter(plccenter==0) = 1;
        
        plccenterref = cclustref(:,cclustID.plcfld);
        plccenterref = round(plccenterref);
        plccenterref(plccenterref==0) = 1;
             
        plccenter = (plccenter-plccenterref);
        plccenter(plccenter<-750) = plccenter(plccenter<-750)+1500;   
        plccenter(plccenter>750) = plccenter(plccenter>750)-1500;
        plccenter = abs(plccenter);
        
        stable =  isplace & isplaceref & plccenter <= 250 & isthere;%
        remap = false(size(stable));%plccenter > 250 & isplace & isplaceref & isthere;
        denovo = isplace & ~isplaceref & isthere;
        disap = ~isplace & isplaceref & isthere;
       
        if j<5
            Nv(:,i) = Nv(:,i) + [sum(stable);sum(denovo);sum(isthere)];
        else
            Nd(:,i) = Nd(:,i) + [sum(stable);sum(denovo);sum(isthere)];
        end
        %% read out binned speed values
        
        speedbin = CAIMcorr(numexp(i+1),j).behave.speedbin;
        speedbin = speedbin(:,1:30);%*100;
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
        
        %% read out place preference vectors
        
        plcvctang = cclust(:,cclustID.plcvctang);
        plcvctang = -(plcvctang-pi);
        plcvct = cclust(:,cclustID.plcvct);

        d = [];
        for k = 1:size(plccenter,1)   % loop over cells        
            d = [d;plcvct(k)*exp(plcvctang(k)*1j)];
        end
        
        if j<5
            PVvStbl = [PVvStbl; d(stable)];
            PVvRemap = [PVvRemap; d(remap)];
            PVvDeno = [PVvDeno; d(denovo)];
        else
            PVdStbl = [PVdStbl; d(stable)];
            PVdRemap = [PVdRemap; d(remap)];
            PVdDeno = [PVdDeno; d(denovo)];
        end
        
        
        %%
        plcfield = CAIMcorr(numexp(i+1),j).plcfield;
        plcfield = plcfield(:,1:150);
        
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
        
        plcstable = plcfield(stable,:);       
        plccenter = cclust(:,cclustID.plcfld);
        [~,h] = sort(plccenter(stable));
        plcstable = plcstable(h,:);
        
        plcremap = plcfield(remap,:);       
        plccenter = cclust(:,cclustID.plcfld);
        [~,h] = sort(plccenter(remap));
        plcremap = plcremap(h,:);
        
        plcdenovo = plcfield(denovo,:);       
        plccenter = cclust(:,cclustID.plcfld);
        [~,h] = sort(plccenter(denovo));
        plcdenovo = plcdenovo(h,:);
        
        if j<5
            PFvStbl = [PFvStbl; plcstable];
            PFvRemap = [PFvRemap; plcremap];
            PFvDeno = [PFvDeno; plcdenovo];
        else
            PFdStbl = [PFdStbl; plcstable];
            PFdRemap = [PFdRemap; plcremap];
            PFdDeno = [PFdDeno; plcdenovo];
        end
        %% correlation coefficients
        r = [0 0];
        if ~isempty(plcstable)
            r(1) = corr(speedbin(1,:)',mean(plcstable,1)');
        end        
%         if ~isempty(plcremap)
%             r(2) = corr(speedbin(1,:)',mean(plcremap,1)');
%         end
        if ~isempty(plcdenovo)
            r(2) = corr(speedbin(1,:)',mean(plcdenovo,1)');
        end      
        r = round(r,2);
        
        if j<5
            R(j,1,:) = r;
        else
            R(j-4,2,:) = r;
        end

        %%
        if DoPlot(2) > 0
            ind = find(stable);
            
            if ~isempty(ind)
                if DoPlot(2) > 1
                    endPlot = length(ind);
                else
                    endPlot = 1;
                end
                for k = 1:endPlot
                    figure('color',[1 1 1],...
                        'renderer','painters',...
                        'visible','on',...
                        'Units','centimeters',...
                        'position',[3 5 [15 15]],...
                        'PaperUnits','centimeters',...
                        'PaperSize', [15 15])
                    for kk = 1:2
                        subplot(5,2,kk:2:7+kk)
                        axon = CAIMcorr(i+kk-1,j).B(:,ind(k));
                        axon = reshape(full(axon),size(CAIMcorr(i+kk-1,j).Cn));
                        imagesc(axon)
                        colormap(jet)
                        title('stable')
                        axis off
                        subplot(5,2,8+kk)
                        imagesc(CAIMcorr(i+kk-1,j).plcfield(ind(k),:))
                        colormap(jet)
                        title('field')
                    end
                end
            end
            ind = find(denovo);
            if ~isempty(ind)
                if DoPlot(2) > 1
                    endPlot = length(ind);
                else
                    endPlot = 1;
                end
                for k = 1:endPlot
                    figure('color',[1 1 1],...
                        'renderer','painters',...
                        'visible','on',...
                        'Units','centimeters',...
                        'position',[3 5 [15 15]],...
                        'PaperUnits','centimeters',...
                        'PaperSize', [15 15])
                    for kk = 1:2
                        subplot(5,2,kk:2:7+kk)
                        axon = CAIMcorr(i+kk-1,j).B(:,ind(k));
                        axon = reshape(full(axon),size(CAIMcorr(i+kk-1,j).Cn));
                        imagesc(axon)
                        title('newly coding')
                        colormap(jet)
                        axis off
                        subplot(5,2,8+kk)
                        imagesc(CAIMcorr(i+kk-1,j).plcfield(ind(k),:))
                        colormap(jet)
                        title('field')
                    end
                end 
            end                                                                               
        end
    end
end

%% Mean Polar plots
if DoPlot(1) == 1
    % % rayleigh test
    % pR(1,1) = circ_rtest(pi+angle(PVvStbl));
    % pR(1,2) = circ_rtest(pi+angle(PVvDeno));
    % pR(2,1) = circ_rtest(pi+angle(PVdStbl));
    % pR(2,2) = circ_rtest(pi+angle(PVdDeno));
    % pR = round(pR,3);    
    % % Herman rasson test
    % pHR(1,1) = hrtest(pi+angle(PVvStbl));
    % pHR(1,2) = hrtest(pi+angle(PVvDeno));
    % pHR(2,1) = hrtest(pi+angle(PVdStbl));
    % pHR(2,2) = hrtest(pi+angle(PVdDeno));
    % pHR = round(pHR,3);    
        
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[3 5 [15 15]],...
        'PaperUnits','centimeters',...
        'PaperSize', [15 15])
    ftsz = 6;
    
    subplot(2,2,1)
    for j = 1:length(PVvStbl)
        polarplot([0 real(-1j*log(PVvStbl(j)))],[0 abs(PVvStbl(j))],'color',[.5 .5 .5]);
        hold on
    end
    d = sum(PVvStbl)/length(PVvStbl);
    
    polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
    polarplot([0 real(-1j*log(Z2))],[0 abs(Z2)],'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(Z3))],[0 abs(Z3)],'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(d))],[0 abs(d)],'color',[1 0 0],"linewidth",2);
    grid off
    ax = gca;
    ax.FontSize = ftsz;
    ax.RTickLabel = {};
    % ax.ThetaTick = [0 60 120 240];
    % ax.ThetaTickLabel = {'Zone2' 'AP' 'Zone1' 'Zone3' };
    title(['I-D, stable coding']);
    % title(['Stable Ventral p_{HR} = ' num2str(pHR(1,1)) ', p_{Ra} = ' num2str(pR(1,1))]);
    
    subplot(2,2,2)
    for j = 1:length(PVvDeno)
        polarplot([0 real(-1j*log(PVvDeno(j)))],[0 abs(PVvDeno(j))],'color',[.5 .5 .5]);
        hold on
    end
    d = sum(PVvDeno)/length(PVvDeno);
    polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
    polarplot([0 real(-1j*log(Z2))],[0 abs(Z2)],'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(Z3))],[0 abs(Z3)],'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(d))],[0 abs(d)],'color',[1 0 0],"linewidth",2);
    grid off
    ax = gca;
    ax.FontSize = ftsz;
    ax.RTickLabel = {};
    ax.ThetaTick = [0 60 120 240];
    ax.ThetaTickLabel = {'Zone2' 'AP' 'Zone1' 'Zone3' };
    title(['I-D, newly coding']);
    % title(['I-D, newly coding, p_{HR} = ' num2str(pHR(1,2)) ', p_{Ra} = ' num2str(pR(1,2))]);
    
    subplot(2,2,3)
    for j = 1:length(PVdStbl)
        polarplot([0 real(-1j*log(PVdStbl(j)))],[0 abs(PVdStbl(j))],'color',[.5 .5 .5]);
        hold on
    end
    d = sum(PVdStbl)/length(PVdStbl);
    polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
    polarplot([0 real(-1j*log(Z2))],[0 abs(Z2)],'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(Z3))],[0 abs(Z3)],'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(d))],[0 abs(d)],'color',[1 0 0],"linewidth",2);
    grid off
    ax = gca;
    ax.FontSize = ftsz;
    ax.RTickLabel = {};
    ax.ThetaTick = [0 60 120 240];
    ax.ThetaTickLabel = {'Zone2' 'AP' 'Zone1' 'Zone3' };
    title(['D-D, newly coding']);
    % title(['D-D, newly coding, p_{HR} = ' num2str(pHR(2,1)) ', p_{Ra} = ' num2str(pR(2,1))]);
    
    subplot(2,2,4)
    for j = 1:length(PVdDeno)
        polarplot([0 real(-1j*log(PVdDeno(j)))],[0 abs(PVdDeno(j))],'color',[.5 .5 .5]);
        hold on
    end
    d = sum(PVdDeno)/length(PVdDeno);
    polarplot([0 real(-1j*log(DZ))],[0 abs(DZ)],'color',[.8 .1 .1]);
    polarplot([0 real(-1j*log(Z2))],[0 abs(Z2)],'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(Z3))],[0 abs(Z3)],'color',[.1 .1 .8]);
    polarplot([0 real(-1j*log(d))],[0 abs(d)],'color',[1 0 0],"linewidth",2);
    grid off
    ax = gca;
    ax.FontSize = ftsz;
    ax.RTickLabel = {};
    ax.ThetaTick = [0 60 120 240];
    ax.ThetaTickLabel = {'Zone2' 'AP' 'Zone1' 'Zone3' };
    title(['D-D, stable coding']);
    % title(['D-D, stable coding p_{HR} = ' num2str(pHR(2,2)) ', p_{Ra} = ' num2str(pR(2,2))]);
end
%% Place preference Vector properties of tracked axons
% length
PlacePrefVect(1,:,i) = [abs(mean(PVvStbl)) abs(mean(PVvDeno)) abs(mean(PVdStbl)) abs(mean(PVdDeno))];
% Std of length
PlacePrefVect(2,:,i) = [std(abs(PVvStbl)) std(abs(PVvDeno)) std(abs(PVdStbl)) std(abs(PVdDeno))];
% SEM of length
PlacePrefVect(3,:,i) = PlacePrefVect(2,:,i)./sqrt([length(PVvStbl) length(PVvDeno) length(PVdStbl) length(PVdDeno)]);
% Angle difference to AP
dist = [(mean(PVvStbl)) (mean(PVvDeno)) (mean(PVdStbl)) (mean(PVdDeno))];
% shift = angle(DZ);%pi;
shift = pi;
shift = exp(shift*1j);
dist = dist./shift;
dist = angle(dist);
dist = -dist;
dist = dist/pi*180;
PlacePrefVect(4,:,i) = dist;
% STD of Angles
PlacePrefVect(5,:,i) = [std(PVvStbl) std(PVvDeno) std(PVdStbl) std(PVdDeno)]/pi*180;
PlacePrefVect(6,:,i) = PlacePrefVect(5,:,i)./sqrt([length(PVvStbl) length(PVvDeno) length(PVdStbl) length(PVdDeno)]);
% Fractions of feature coders
PlacePrefVect(7,:,i) = [Nv(1,i)/Nv(3,i) Nv(2,i)/Nv(3,i) Nd(1,i)/Nd(3,i) Nd(2,i)/Nd(3,i)];


end
