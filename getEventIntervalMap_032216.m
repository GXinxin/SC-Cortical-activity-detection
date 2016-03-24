% compute roi event interval for SC and/or V1 recordings
% need roi contours, and several movies

% load('15061301_regress.mat');
% clear; clc;

% cd 'D:\Lab\Data\ChAT\160115'
% cd 'D:\Lab\Data\ChAT\160317'


addpath /fastscratch/xg77/CHAT/160317
addpath /home2/xg77/CalciumImgCode  
addpath(genpath('/home2/xg77/CalciumImgCode/piotr_toolbox'))
addpath(genpath('/home2/xg77/CalciumImgCode/wholeBrainDX'))
addpath(genpath('/home2/xg77/CalciumImgCode/sigTOOL'))
addpath(genpath('/home2/xg77/CalciumImgCode/CalciumDX'))
addpath(genpath('/home2/xg77/CalciumImgCode/bfmatlab'))
addpath(genpath('/home2/xg77/CalciumImgCode/chatAnalysis'))
addpath(genpath('/home2/xg77/CalciumImgCode/piotr_toolbox'))



filelist = readtext('files.txt', ' ');
fnms = filelist(:, 1);
concat_fnm = filelist(:, 2);
mask_fnms = filelist(:, 4);
svd_fnms = filelist(:, 3);
is_stim = filelist(:, 5);
no_movies = length(fnms);
downSampleRatio = 0.5;
sigma = 3; % sigma for gaussian smoothing
thresh = 5; % segmentation threshold

for n = 1:no_movies
    
    % load masks for SC and V1
    load(mask_fnms{n});
    
    % load movie (concat two movies together, 6000 frames in total)
    fnm = fnms{n};
    img1 = openMovie(fnm);
    img1 = imresize(img1, downSampleRatio, 'bilinear');
    img2 = openMovie(concat_fnm{n});
    img2 = imresize(img2, downSampleRatio, 'bilinear');
    imgall = cat(3, img1, img2);
       
    
    % compute dF/F
    sz = size(imgall);
    imgall = reshape(imgall, sz(1) * sz(2), sz(3));
    
    Amean = mean(imgall, 2);
    imgall = imgall ./ (Amean * ones(1, size(imgall, 2))) - 1;
    
    
    % apply SVD, reconstruct based on good PCs
    load(svd_fnms{n});
    eigenLoad = imgall' * mixedfilters2;
    imgall = eigenLoad * mixedfilters2';
    imgall = imgall';
    imgall = zscore(reshape(imgall, 1, sz(1)*sz(2)*sz(3)));
    
    
    % apply gaussian smooth
    imgall = reshape(imgall, sz(1), sz(2), sz(3));
    parfor fr = 1:sz(3); 
        imgall(:, :, fr) = gaussSmooth(imgall(:,:,fr), sigma, 'same');
    end	
     
   
    totalMask = roi{1} + roi{2} + roi{3};
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute flow field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imgall = reshape(imgall, sz(1), sz(2), sz(3));   
    [normVx, normVy] = computeFlowField_normalized_xx(imgall, sz, 0.25);
    
    smallMask = imresize(totalMask, 0.125, 'bilinear');
    h = figure; imagesc(smallMask); hold on
    quiver(mean(normVx, 3).*smallMask, mean(normVy, 3).*smallMask);
    title(fnm(1:end-4))
%     set(h, 'Position', [0, 0, 2000, 2000]);
%     h.PaperPositionMode = 'auto';
%     print('ScreenSizeFigure','-dpng','-r0')
%     set(h,'PaperPositionMode','auto')
    print(h, [fnm(1:end-4), '_quiver'], '-dpdf')
%     print('-bestfit','BestFitFigure','-dpdf')
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(roi)
%     for r = 1:2
        
        savefn = [fnm(1:end-4), '_', num2str(r)];
        
        maskMatrix = [];
        if r == length(roi) + 1
            maskMatrix = imresize(totalMask, downSampleRatio, 'bilinear');
        else
            maskMatrix = imresize(roi{r}, downSampleRatio, 'bilinear');
        end
        maskId = find(maskMatrix > 0);
        
        maskMatrix = reshape(maskMatrix, sz(1)*sz(2), 1);
        maskMatrix = repmat(maskMatrix, 1, sz(3));
%         save('check.mat', 'maskMatrix', 'sz')
        maskMatrix = reshape(maskMatrix, sz(1), sz(2), sz(3));        
        

        imgall = reshape(imgall, sz(1), sz(2), sz(3));
        subMov = imgall .* maskMatrix;

        
        % zscore
        subMov = zscore(reshape(subMov, 1, sz(1)*sz(2)*sz(3)));
        if r == length(roi)
            thresh = 4;
        end
        activeMov = subMov > thresh;
        subMov = reshape(subMov, sz(1), sz(2), sz(3));
        activeMov = reshape(activeMov,sz(1), sz(2), sz(3));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % detect connected components
        CC = bwconncomp(activeMov);
        STATS = regionprops(CC, activeMov);  %some of the properties in regionprops that work on n-D arrays


        % compute duration, diameter, area and centroid
        roiBoundingBox = zeros(length(STATS),6);
        for i = 1:length(STATS)
            roiBoundingBox(i, :) = STATS(i).BoundingBox;
        end

        dura = roiBoundingBox(:,6); 
        dia = mean([roiBoundingBox(:,4) roiBoundingBox(:,5)], 2);
        area = STATS.Area;
        
        if r < length(roi)
            validId{r, n} = (dura > 30 & dia > 30); % this works for cAMP SC data
        else
            validId{r, n} = (dura > 2 & dia > 20 & area > 200); % this works ok for cortex, 160317_rightSpon_cAMP31
        end
        durations{r, n} = dura(validId{r, n});
        diameters{r, n} = dia(validId{r, n});
        roiCentr{r, n} = vertcat(STATS(validId{r, n}).Centroid);
        roiArea{r, n} = [STATS(validId{r, n}).Area];  

        
        % Compute wave direction
        [AVx, AVy] = computeFlowField_xx(imgall, sz); 
              
        nDomains = sum(validId{r, n});
        Vsum{r}(1:nDomains) = struct('theta',[],'rho',[]);
        validDomains = CC.PixelIdxList(validId{r, n});
        
        for k = 1:nDomains
            [theta, rho]= cart2pol(sum(AVx(validDomains{k})), sum(AVy(validDomains{k})));
            Vsum{r}(k).theta = theta;
            Vsum{r}(k).rho = rho;
        end

        angle{r}{n} = vertcat(Vsum{r}.theta);
        h = figure; rose(angle{r}{n});
        set(gca,'YDir','reverse'); % 90 degree is moving downwards
        title(savefn);
        saveas(h, [savefn, '_rosePlot.png'])
        

        % plot durations and diameters of detected components
        % savefn = ['dA_thresh', num2str(thresh1), 'sigma_', num2str(sigma)];            
%         h(1) = figure; hist(durations{n}, 50); xlabel('durations (frames)'); title(['thresh=', num2str(thresh)])
        % saveas(h(1), [savefn, '_durations.png'])
%         h(2) = figure; hist(diameters{n}, 50); xlabel('diameters (pixels)'); title(['thresh=', num2str(thresh)])
        % saveas(h(2), [savefn, '_diameters.png'])
%         h(3) = figure; scatter(durations{n}, diameters{n}); xlabel('durations'); ylabel('diameters'); title(['thresh=', num2str(thresh)])
        % saveas(h(3), [savefn, '_duraVSdia.png'])  


        % reconstruct binary mov based on valid connected components
        valid_mov = subMov;
        valid_activeMov = activeMov;
        badId = find(validId{r, n} == 0);
        for id = 1:length(badId)
            removePixel = CC.PixelIdxList{badId(id)};
            valid_mov(removePixel) = 0;
            valid_activeMov(removePixel) = 0;
        end

        valid_activeMov_down = imresize(valid_activeMov, 0.25, 'bilinear');
        tmpsz = size(valid_activeMov_down);
        valid_activeMov_down = reshape(valid_activeMov_down, tmpsz(1)*tmpsz(2), tmpsz(3));



        % create segmented movie
        Iarr=mat2gray(imgall);   %scale the whole array
%         [I2arr, map] = gray2ind(Iarr, 256); %convert the whole array to 8bit indexed		
        for fr = 1:sz(3)
%             M(fr) = im2frame(I2arr(:,:,fr), map);  %setup the indexed raw dFoF movie
            [I2, map2] = gray2ind(valid_activeMov(:,:,fr), 8); %figure; imshow(I2,map)
            F(fr) = im2frame(I2,map2);  %setup the binary segmented mask movie
        end
    %     fnm2 = [savefn, '_smoothed_dFoF.avi']; 
    %     writeMovie(M,fnm2);
        fnm3 = [savefn, '_mask.avi'];
        writeMovie(F,fnm3);


        
        
        % compute active duration and event interval   
        for p = 1:tmpsz(1)*tmpsz(2)
            activeOn{p} = find(valid_activeMov_down(p, 2:end) - valid_activeMov_down(p, 1:end-1) > 0) + 1;
            activeOff{p} = find(valid_activeMov_down(p, 2:end) - valid_activeMov_down(p, 1:end-1) < 0);

            if (isempty(activeOn{p}) + isempty(activeOff{p})) == 1

                activeOn{p} = [];
                activeOff{p} = [];

            elseif (isempty(activeOn{p}) + isempty(activeOff{p})) == 0

                if activeOn{p}(1) > activeOff{p}(1)
                    activeOff{p} = activeOff{p}(2:end);
                end

                if activeOn{p}(end) > activeOff{p}(end)
                    activeOn{p} = activeOn{p}(1:end-1);
                end
            end

            pixelDuration{p} = activeOff{p} - activeOn{p};

            goodId = pixelDuration{p} > 2;
            activeOn{p} = activeOn{p}(goodId);
            activeOff{p} = activeOff{p}(goodId);
            pixelDuration{p} = pixelDuration{p}(goodId);

            pixelInterval{p} = activeOn{p}(2:end) - activeOff{p}(1:end-1);
            meanInterval(p) = mean(pixelInterval{p});
        end
        
        p_Interval{r, n} = pixelInterval;
        
        meanInterval(isnan(meanInterval)) = 2000;

        meanInterval = reshape(meanInterval, tmpsz(1), tmpsz(2));
        h = figure; imagesc(meanInterval); colorbar; caxis([0, 2200]); colormap jet
        title(savefn);
        m_Interval{r, n} = meanInterval;
        
        saveas(h, [savefn, '_interval.png']);
        
    end
    
    
end
save('dataSummary_rightSpon.mat', 'fnms', 'p_Interval', 'm_Interval', 'validId', 'durations', 'diameters', 'roiCentr', 'roiArea', 'angle');


