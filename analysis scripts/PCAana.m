function PCAout = PCAana(caim,scn,method)
if ~isfield(caim,'S')
    PCAout = [];
    return
end
addpath(genpath('/media/2Photon/Matlab code/Martin Source Code/gpfa_v0203'))
% C = caim.C;C = C./caim.Df;  

C = caim.S;
% C(caim.S_bin==0)= 0;

for i = 1:size(C,1)
%     C(i,:) = zscore(C(i,:));
%     C(i,:) = C(i,:) + abs(min(C(i,:)));
    C(i,:) = C(i,:) /max(C(i,:) );
end

clear dat
dat.trialId = 1;
dat.spikes = C;
%%
% [a,score,~,~,explained] = pca(C);
%%
% Select method to extract neural trajectories:
% 'gpfa' -- Gaussian-process factor analysis
% 'fa'   -- Smooth and factor analysis
% 'ppca' -- Smooth and probabilistic principal components analysis
% 'pca'  -- Smooth and principal components analysis
if strcmp(method,'pca')
    disp('Good old PCA')
    runIdx = 2;
    method = 'pca';
    xDim = size(C,1)-1;
    binWidth = 1; 
    kernSD = 10;

    % Extract neural trajectories
    result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                        'kernSDList', kernSD,'binWidth',binWidth);
    if isempty(result);PCAout = [];return;end
    [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);%%
    score = result.kern.estParams.L;
    hasSpikesBool = result.hasSpikesBool;
    if sum(hasSpikesBool)<xDim
        xDim = sum(hasSpikesBool);
    end
%    score = estParams.Corth;
    a = [];
    for i = 1:size(seqTrain,2)
        a = [a seqTrain(i).xorth];
    end
    clear b

    for i = 1:size(a,1)
        b(i,:) = interp1(1:size(a,2),a(i,:),1:size(a,2)/size(caim.C,2):size(a,2),'spline');
    end

    b(:,end+1:size(caim.C,2)) = 0;
    a = b';
elseif strcmp(method,'ica')
    %%
    disp('I am performing a really nice ICA, bitches')
    kernSD = 5;
    xDim = 30;
    yOut = smoother(C, kernSD,1);
    b = rica(yOut',xDim,'IterationLimit',1000);
    % b = rica(a,q,'IterationLimit',1000);
    a = yOut'*b.TransformWeights;
    score = b.TransformWeights;
    hasSpikesBool = true(size(C,1),1);
end


%%
explained = var(a);

linEqn = 'a+b*(x)';
int = 1:round(length(explained)/3);
% x = log(1:length(int));
x = log(int);
y = log(explained(int));
startPoints = [y(1) -1 ];
upperBounds = [y(1)+3 1];
lowerBounds = [y(1)-3 -3 ];
f1 = fit(x',y',linEqn,...
       'Start', startPoints,...
       'Upper',upperBounds,...
       'Lower',lowerBounds);
   Alpha = -round(f1.b,2);
int = 1:length(explained);
x = log(1:length(int));
y = log(explained(int));

PCAout.alpha = Alpha;
PCAout.explained = explained;


%% Sectioned PCA
% if strcmp(method,'pca') && sum(scn.running)>150 
%     PCAout.sectioned = PCAsectioned(caim,scn,'pca');    
% else
%     PCAout.sectioned = [];
% end

end