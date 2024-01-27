function TrackedFOVPlot(CAIMcorr,cclustID)
    mouseID = {'M259','M261', 'M270','M272','M262','M263','M271','M278'};
    
    numPlt = size(CAIMcorr,1);
    numMic = size(CAIMcorr,2);
    
    % scaling factors for the individual images to adjust brightness
    scl = [2 2 5 6 2 4 4 5;
            2 2 5 4 2 4 3 5;
             2 1 4 4 2 4 3 5];
         
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[3 5 [15 30]],...
        'PaperUnits','centimeters',...
        'PaperSize', [15 30])
    
    for j = 1:numMic
        A = CAIMcorr(1, j).Cncor;
        
        for i = 1:numPlt
            A(:,:,i) = CAIMcorr(i,j).Cncor;
            subplot(numMic, numPlt+1,(j-1)*(numPlt+1)+i)
            imshow(scl(i,j)*mat2gray(A(:,:,i)))
            if i == 1
                title([mouseID{j} ', Day1'])
            else
                title(['Day2'])
            end
        end
    
        C = imfuse(A(:,:,1),A(:,:,2));
        subplot(numMic,numPlt+1,j*(numPlt+1))
        imshow(scl(i+1,j)*C)
        title(['Composite'])
    end
    
    %%
    
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[3 5 [25 60]],...
        'PaperUnits','centimeters',...
        'PaperSize', [25 60])
    
    % scaling factors for the individual images to adjust brightness
    % Green channel
    scl = [14 1 5 6 6 4 4 3;
            8 1 3 7 6 4 7 7];
    
    % red channel
    scl(:,:,2) = [3 2 7 8 5 5 3 7;
                   4 3 7 7 5 5 3 8];
    for j = 1:numMic
        for i = 1:numPlt
            FOV = scl(i,j,2)*mat2gray(CAIMcorr(i,j).FOV(:,:,2));
            FOV(:,:,2) = scl(i,j,1)*mat2gray(CAIMcorr(i,j).FOV(:,:,1));
            FOV(:,:,3) = 0;
            
            if i >1
                sameAsInput = affineOutputView(size(FOV(:,:,1)),CAIMcorr(i,j).tform,'BoundsStyle','SameAsInput');
                FOV(:,:,1) = imwarp(FOV(:,:,1),CAIMcorr(i,j).tform,'OutputView',sameAsInput);
                FOV(:,:,2) = imwarp(FOV(:,:,2),CAIMcorr(i,j).tform,'OutputView',sameAsInput);
            end
            subplot(numMic, numPlt,(j-1)*(numPlt)+i)
            imshow(FOV)
            if i == 1
                title([mouseID{j} ', Day1'])
            else
                title(['Day2'])
            end
        end
        
    end

end