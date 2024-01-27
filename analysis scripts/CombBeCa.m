function [CAIM,PCA,caim,scn,behave,TRACE] = CombBeCa(pathname,mouseID)

[files,samedate,prefiles] = FindDataSets(pathname);

% behave.expID = cell2mat(samedate(:,2));
behave.expID = 1:length(files);
behave.numrounds = zeros(length(files),1);
behave.runtime = zeros(length(files),1);
behave.pupilsize = [];
behave.distance = []; %Negar added this 

CAIM = [];
PCA = [];
TRACE = [];
%%
for i = 1:length(files)
    close all
    %% Load Data if already processed before
    if ~isempty(prefiles) && length(prefiles)>=i
        disp(['load ' prefiles(i).name])
        load(prefiles(i).name,'belt','caim','scn');  
    else
        %% Readout ca-Imaging & behaviour data and real time correction  
        [belt,caim] = readcaim(pathname,files(i).name(1:end-6));
        % !!This function concatenates experiements that happened on the
        % same day!!
%         [belt,caim] = ConnecCaim(pathname,files,samedate{i,1});   
        [belt,scn] = BeltToSCN(caim,belt);   
    end
    if isempty(caim.Cn)
        Cn = zeros(caim.options.d1,caim.options.d2,size(caim.A,2));
        for ii = 1:size(caim.A,2)
            Cn(:,:,ii) = full(reshape(caim.A(:,ii),caim.options.d1,caim.options.d2));
        end
        caim.Cn = mean(Cn,3);
    end
    %% Place and Speed coding analysis
     scn = placecor(caim,scn);
    %% Firing properties analysis
     caim.fireprop = firefreq(caim,scn);  %comment in & out
    
    %% Ensemble analysis  
    runshuf = 1; % comment in for the first time you analyse something
%     comment out when you want to run the whole structure after the first
%     time 
    % Similarity measures logicals: shuf, pca, binary, kernSD, norm  

%     caim.PCAout = PCAana(caim,scn,'pca'); %comment in & out
    caim.GLMout = GLMana(caim,scn,'pca', runshuf);
%     caim.ICAout = PCAana(caim,scn,'ica');
%     caim.GPFAout = PCAsectioned(caim,scn,'gpfa');
    %% Input Bulk Signal analysis
%     caim = bulkanalysis(caim,scn);  %comment in & out    
    %% Correlation to stimuli
%      scn = stimcor(belt,caim,scn);  %comment in & out   
    %% Pool Data
    %  Behaviour
    behave.numrounds(i) = max(scn.rounds);
    behave.runtime(i) = scn.tsscn(end)*sum(scn.running==1)/length(scn.tsscn)/1000;
    behave.pupilsize = cat(3,behave.pupilsize,scn.pupilsize);  
   
    %% Cell properties   
    scn.cclust = cellcluster(caim,scn,str2double(mouseID(2:end)));     

    %% Pooling struct output
    
    if isfield(scn,'cellID') && size(scn.spaccode,2)>1
        % Pool placefield plot (3D variable. slice 1 is raw data, slice 5
        % is the result after Dombeck arguments)
        plcfield = NaN(size(scn.spaccode,1),size(scn.spaccode,2));
        plcfield(scn.cellID(:,1),:) = scn.plcfield(:,:,1);
        overround = scn.overround;
        restcode = mean(scn.restcode,3);
    elseif isfield(caim,'Y')
        plcfield = NaN(size(caim.Y,1),150);
        overround = NaN(size(caim.Y,1),1);
    else
        plcfield = [];
        overround = [];
    end
    
    CAIM(behave.expID(i)).behave.numrounds = max(scn.rounds);
    CAIM(behave.expID(i)).behave.totdist = scn.totdist(end-1);
    CAIM(behave.expID(i)).behave.meanspeed = 100*mean(scn.speed(scn.running==1));
    CAIM(behave.expID(i)).behave.runtime = scn.tsscn(end)*sum(scn.running==1)/length(scn.tsscn)/1000;
    CAIM(behave.expID(i)).behave.pupilsize = scn.pupilsize;
    CAIM(behave.expID(i)).behave.running = logical(scn.running);
    CAIM(behave.expID(i)).behave.distance = scn.distance;
    CAIM(behave.expID(i)).behave.speed = scn.speed;
    CAIM(behave.expID(i)).behave.tsscn = scn.tsscn;
    CAIM(behave.expID(i)).behave.speedwaitbin = scn.speedwaitbin;
    CAIM(behave.expID(i)).behave.speedbin = scn.speedbin;
    CAIM(behave.expID(i)).behave.waitbin = scn.waitbin;
    CAIM(behave.expID(i)).behave.distance = scn.distance;
    if isfield(caim,'Y')
        CAIM(behave.expID(i)).Y = [];
        CAIM(behave.expID(i)).A = caim.A;
        CAIM(behave.expID(i)).S = logical(caim.S_bin);
        CAIM(behave.expID(i)).SRaw = single(caim.S);
        CAIM(behave.expID(i)).C = single(caim.C);
        CAIM(behave.expID(i)).Cn = caim.Cn; 
        CAIM(behave.expID(i)).FOV = caim.FOV;
        CAIM(behave.expID(i)).cclust = scn.cclust;
        CAIM(behave.expID(i)).plcfield = plcfield;
        CAIM(behave.expID(i)).overround = overround;   
        CAIM(behave.expID(i)).restcode = restcode;   
%         CAIM(behave.expID(i)).network = caim.network;
        CAIM(behave.expID(i)).speedcorr = scn.speedcorr;
        CAIM(behave.expID(i)).fireprop = caim.fireprop;
        CAIM(behave.expID(i)).stripe1 = scn.stripe1;
        CAIM(behave.expID(i)).stripe2 = scn.stripe2;
        CAIM(behave.expID(i)).stripe3 = scn.stripe3;
        PCA(behave.expID(i)).PCAout = caim.PCAout;
%         PCA(behave.expID(i)).ICAout = caim.ICAout;
%         PCA(behave.expID(i)).GPFAout = caim.GPFAout;
        PCA(behave.expID(i)).GLMout = caim.GLMout;
        if isfield(scn,'airpuff')
            CAIM(behave.expID(i)).AP = scn.airpuff.stimon;
        else
            CAIM(behave.expID(i)).AP = false(size(scn.tsscn));
        end
        TRACE(behave.expID(i)).S = single(caim.S);%negar added this; Martin changed it
        TRACE(behave.expID(i)).C = single(caim.C);
    else
        
        CAIM(behave.expID(i)).Y = [];
        CAIM(behave.expID(i)).A = [];
        CAIM(behave.expID(i)).S = [];
        CAIM(behave.expID(i)).Cn = [];
        CAIM(behave.expID(i)).FOV = [];
        CAIM(behave.expID(i)).cclust = [];
        CAIM(behave.expID(i)).plcfield = [];
        CAIM(behave.expID(i)).overround = [];
        CAIM(behave.expID(i)).restcode = [];
%         CAIM(behave.expID(i)).network = [];
        CAIM(behave.expID(i)).speedcorr = [];
        CAIM(behave.expID(i)).fireprop = [];
        CAIM(behave.expID(i)).stripe1 = [];
        CAIM(behave.expID(i)).stripe2 = [];
        CAIM(behave.expID(i)).stripe3 = [];
        PCA(behave.expID(i)).PCAout = [];
%         PCA(behave.expID(i)).ICAout = [];
%         PCA(behave.expID(i)).GPFAout = [];
        PCA(behave.expID(i)).GLMout = [];
        TRACE(behave.expID(i)).S = [];
        TRACE(behave.expID(i)).C = [];
    end
    
    if isfield(caim,'bulk')
        CAIM(behave.expID(i)).bulk = caim.bulk;
    else
        CAIM(behave.expID(i)).bulk = [];
    end   
    %% save preproceesed data
    if ispc
        if ~isdir([pathname '\preprocessed']);mkdir([pathname '\preprocessed']);end
        save([pathname '\preprocessed\' files(i).name(1:end-6) 'pro.mat'],'caim','belt','scn')
    else
        if ~isdir([pathname '/preprocessed']);mkdir([pathname '/preprocessed']);end
        save([pathname 'preprocessed/' files(i).name(1:end-6) 'pro.mat'],'caim','belt','scn')
    end
    %% Summary Plot and .pdf
%     SummaryPlot(pathname,files(samedate{i,1}(1)).name(1:end-6),caim,belt,scn)
end

end



