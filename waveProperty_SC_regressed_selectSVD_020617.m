% compute roi event interval for SC and/or V1 recordings
% corrected angle from 032216 version
% add mean active duration for each voxel
% need roi contours, and several movies

% changed interval definition: 
% (previous) on - off
% (new) inter-eventCenter-interval
% 06/22/16


clear; clc;


% cd 'E:\Lab\Data\ChAT\B2\180122_p8_paxB2flox_null_sinusaav9synG6s_4.9g_male'
% cd 'Y:\Swarna_FlailerMice\WT_180216_controlflailer_eyeAAV1G6s_p6_3.1g_110x'


addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/piotr_toolbox'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/wholeBrainDX'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/sigTOOL'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/CalciumDX'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/bfmatlab'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/chatAnalysis'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/piotr_toolbox'))
 


filelist = readtext('files.txt', ' ');
fnms = filelist(:, 1);
mask_fnms = filelist(:, 2);


isSVD = 0;
if isSVD
    svd_fnms = filelist(:, 3);
    thresh = filelist(:, 4);
else
    thresh = filelist(:, 3);
end

no_movies = length(fnms);
downSampleRatio = 1;
maskDownSample = 0.5;
% thresh = 1.7;
sigma = 1; % sigma for gaussian smoothing
dura_th = 8;
dia_th = 8;
sz = [256 250];
% sz = [270 320];
% sz = [251 401];

isTopHat = 1;
diskSz = 20;

frameRate = 10; % Hz


for n = 1:no_movies

    clear ddiff
    
    % load masks for SCs
    load(mask_fnms{n});
    
    % load movie (preprocessed)
    fnm = fnms{n};
    load(fnm)
    if ~exist('ddiff', 'var')
        imgall = dA(:, 50:end);
        clear dA
    else
        imgall = ddiff(:, 50:end);
        clear ddiff
    end
    
    % convert to double if data was saved as single
    if isa(imgall, 'double')
    else
        imgall = double(imgall);
    end
    
    sz(3) = size(imgall, 2);
    imgall(isnan(imgall)) = 0;
    
    adjust = size(imgall, 1) > sz(1)*sz(2);
    if adjust
        imgall = imresize(reshape(imgall, sz(1)*2, sz(2)*2, sz(3)), .5, 'bilinear');
        imgall = reshape(imgall, sz(1)*sz(2), sz(3));
    end
  
    
    % apply SVD, reconstruct based on good PCs
    if isSVD 
        load(svd_fnms{n});
        eigenLoad = imgall' * mixedfilters2;
        imgall = eigenLoad * mixedfilters2';
        imgall = imgall';
    end
    
    % zscore updated to pixel based
%     imgall = zscore(reshape(imgall, 1, sz(1)*sz(2)*sz(3)));
    tmp = zscore(imgall');
    imgall = tmp';
    clear tmp
    
    
    % apply gaussian smooth and/or top-hat filtering
    imgall = reshape(imgall, sz(1), sz(2), sz(3));
    if (isTopHat)
        parfor fr = 1:sz(3) 
            imgall(:, :, fr) = gaussSmooth(imgall(:,:,fr), sigma, 'same');
            se = strel('disk',diskSz);
            imgall(:, :, fr) = imtophat(imgall(:, :, fr), se);
        end	
    else
        parfor fr = 1:sz(3) 
            imgall(:, :, fr) = gaussSmooth(imgall(:,:,fr), sigma, 'same');
        end	
    end

    
    
    totalMask = zeros(size(roi{1}));
    for r = 1:length(roi)
        totalMask = totalMask + roi{r};
    end
    if adjust
        totalMask = imresize(totalMask, .5, 'bilinear');
    end
    totalMask = totalMask > 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute flow field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    smallSize = 0.5; % downSample the flowfield size
    imgall = reshape(imgall, sz(1), sz(2), sz(3));   
    [normVx, normVy] = computeFlowField_normalized_xx(imgall, sz, smallSize);

   
    smallMask = imresize(totalMask, smallSize, 'bilinear');
    
    h = figure; imagesc(smallMask); hold on
    quiver(mean(normVx, 3).*smallMask, mean(normVy, 3).*smallMask); axis image;
    title(fnm(1:end-4))
    set(h, 'Position', [0, 0, 1200, 900]);
    h.PaperPositionMode = 'auto';
    print([fnm(1:end-4), '_quiver'], '-dpng', '-r0')
    
    h = figure; imagesc(imresize(smallMask, .5, 'bilinear')); hold on
    quiver(imresize(mean(normVx, 3), .5, 'bilinear') .* imresize(smallMask, .5, 'bilinear'), ...
        imresize(mean(normVy, 3), .5, 'bilinear') .* imresize(smallMask, .5, 'bilinear'));
    axis image;
    title(fnm(1:end-4))
    set(h, 'Position', [0, 0, 1200, 900]);
    h.PaperPositionMode = 'auto';
    print([fnm(1:end-4), '_quiver_s'], '-dpng', '-r0')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    total_ActiveMovie = [];
    total_AVx = [];
    total_AVy = [];
    for r = 1:length(roi)
        
        savefn = [fnm(1:end-4), '_', num2str(r)];
        
        if r == length(roi) + 1
            maskMatrix = totalMask;
%             maskMatrix = imresize(totalMask, maskDownSample, 'bilinear');
        else
            maskMatrix = roi{r};
%             maskMatrix = imresize(roi{r}, maskDownSample, 'bilinear');
            if adjust
                maskMatrix = imresize(roi{r}, .5, 'bilinear');
            end
    %             mas
        end
        maskId = find(maskMatrix > 0);
        
        maskMatrix = reshape(maskMatrix, sz(1)*sz(2), 1);
        maskMatrix = repmat(maskMatrix, 1, sz(3));
        maskMatrix = reshape(maskMatrix, sz(1), sz(2), sz(3));        
        

        imgall = reshape(imgall, sz(1), sz(2), sz(3));
        subMov = imgall .* maskMatrix;
        
       
        
        % zscore
        subMov = reshape(subMov, sz(1) * sz(2), sz(3));
        temp_subMov = subMov(maskId, :);
        temp_subMov = reshape(zscore(temp_subMov(:)), length(maskId), sz(3));
        subMov = zeros(size(subMov));
        subMov(maskId, :) = temp_subMov;
%         activeMov = subMov > thresh;
        activeMov = subMov > thresh{n};
%         activeMov = (subMov > thresh) + (subMov < -thresh);
        subMov = reshape(subMov, sz(1), sz(2), sz(3));
        activeMov = reshape(activeMov,sz(1), sz(2), sz(3));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % detect connected components
        CC = bwconncomp(activeMov);
        STATS = regionprops(CC, activeMov, 'Area', 'BoundingBox', 'Centroid', 'PixelList');      

        % compute duration, diameter, area and centroid
        roiBoundingBox = zeros(length(STATS),6);
        for i = 1:length(STATS)
            roiBoundingBox(i, :) = STATS(i).BoundingBox;
        end

        dura = roiBoundingBox(:,6); 
        dia = mean([roiBoundingBox(:,4) roiBoundingBox(:,5)], 2);
        area = vertcat(STATS.Area);
        

        validId{r, n} = (dura > dura_th & dia > dia_th); 

        durations{r, n} = dura(validId{r, n});
        diameters{r, n} = dia(validId{r, n});
        roiCentr{r, n} = vertcat(STATS(validId{r, n}).Centroid);
        roiArea{r, n} = [STATS(validId{r, n}).Area];  
        boundBox{r, n} = [STATS(validId{r, n}).BoundingBox];
        valid_tmp = find(validId{r, n} > 0);
        for v = 1:length(valid_tmp)
            pixel{r, n}{v} = STATS(valid_tmp(v)).PixelList;
        end
        
        valid{r, n} = find(validId{r, n} > 0);

        
        
        
        % Compute wave direction
        [AVx, AVy] = computeFlowField_xx(imgall, sz); 
              
        nDomains = sum(validId{r, n});
        validDomains = CC.PixelIdxList(validId{r, n});
        
        theta = []; rho = [];
        for k = 1:nDomains
            p_id = intersect(1 : sz(1)*sz(2)*size(AVx, 3), validDomains{k});
            [theta(k), rho(k)]= cart2pol(sum(AVx(p_id)), sum(AVy(p_id)));
        end
        
        AVx = imresize(AVx, .5, 'bilinear');
        AVy = imresize(AVy, .5, 'bilinear');
        
        if isempty(total_AVx)
            total_AVx = AVx;
            total_AVy = AVy;
        else
            total_AVx = total_AVx + AVx;
            total_AVy = total_AVy + AVy;
        end
        clear AVx AVy
        

        angle{r, n} = theta;
        RHO{r, n} = rho;
        h = figure; rose(angle{r, n});
        set(gca,'YDir','reverse'); % 90 degree is moving downwards
        title(savefn);
        saveas(h, [savefn, '_rosePlot.png'])
        
        

        % plot durations and diameters of detected components
        savefn2 = [savefn, num2str(thresh{n}), 'sigma_', num2str(sigma)];            
        h(1) = figure; hist(durations{r, n}, 50); xlabel('durations (frames)'); title(['thresh=', num2str(thresh{n})])
        saveas(h(1), [savefn2, '_durations.png'])
        h(2) = figure; hist(diameters{r, n}, 50); xlabel('diameters (pixels)'); title(['thresh=', num2str(thresh{n})])
        saveas(h(2), [savefn2, '_diameters.png'])
        h(3) = figure; scatter(durations{r, n}, diameters{r, n}); xlabel('durations'); ylabel('diameters'); title(['thresh=', num2str(thresh{n})])
        saveas(h(3), [savefn2, '_duraVSdia.png'])  


        
        % reconstruct binary mov based on valid connected components
        valid_mov = subMov;
        valid_activeMov{r} = activeMov;
        clear activeMov subMov
        
        badId = find(validId{r, n} == 0);
        for id = 1:length(badId)
            removePixel = CC.PixelIdxList{badId(id)};
            valid_mov(removePixel) = 0;
            valid_activeMov{r}(removePixel) = 0;
        end

        valid_activeMov_down{r} = imresize(valid_activeMov{r}, 0.5, 'bilinear');
        tmpsz = size(valid_activeMov_down{r});
        valid_activeMov_down{r} = reshape(valid_activeMov_down{r}, tmpsz(1)*tmpsz(2), tmpsz(3));
        
        if isempty(total_ActiveMovie)
            total_ActiveMovie = valid_activeMov{r};
        else
            total_ActiveMovie = total_ActiveMovie + valid_activeMov{r};
        end


        
        
        % compute active duration and event interval   
        for p = 1:tmpsz(1)*tmpsz(2)
            activeOn{p} = find(valid_activeMov_down{r}(p, 2:end) - valid_activeMov_down{r}(p, 1:end-1) > 0) + 1;
            activeOff{p} = find(valid_activeMov_down{r}(p, 2:end) - valid_activeMov_down{r}(p, 1:end-1) < 0);

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
            meanDuration(p) = sum(pixelDuration{p}) / sz(3);
            
            % interval: from end of one event to the beginning of the following event
            pixelInterval1{p} = activeOn{p}(2:end) - activeOff{p}(1:end-1);
            meanInterval1(p) = mean(pixelInterval1{p});
            
            % interval: between the center of each event
            eventCenterTime{p} = activeOff{p} - activeOn{p};
            pixelInterval2{p} = eventCenterTime{p}(2 : end) - eventCenterTime{p}(1 : end-1);
            meanInterval2(p) = mean(pixelInterval2{p});
        end

        
        p_Interval1{r, n} = pixelInterval1;       
        meanInterval1(isnan(meanInterval1)) = 2000/frameRate;
        meanInterval1 = reshape(meanInterval1, tmpsz(1), tmpsz(2))/frameRate;
        h = figure; imagesc(meanInterval1); colorbar; colormap jet
        caxis([0, 50]); axis image
        title(savefn2);
        saveas(h, [savefn2, '_interval.png']);       
               
        
        meanDuration = reshape(meanDuration, tmpsz(1), tmpsz(2));
        h = figure; imagesc(meanDuration); colorbar; colormap jet
        caxis([0, 0.3]); axis image
        title(savefn2);
        saveas(h, [savefn2, '_duration.png']);
        
        m_Interval1{r, n} = meanInterval1;
        m_p_Duration{r, n} = meanDuration;
                     
    end 
    clear valid_activeMov imgall
    
    
    % create segmented movie     
    total_ActiveMovie = total_ActiveMovie > 0;
    for fr = 1:sz(3)
        [I2, map2] = gray2ind(total_ActiveMovie(:,:,fr), 8); %figure; imshow(I2,map)
        F(fr) = im2frame(I2,map2);  %setup the binary segmented mask movie
    end
    fnm3 = [fnm(1:end-4), '_mask_th', num2str(thresh{n}), '.avi'];
    writeMovie_xx(F, fnm3, 0);
    
    
    total_ActiveMovie = imresize(total_ActiveMovie, .5, 'bilinear');
    
    save([fnm(1:end-4), '_opticFlow.mat'], 'fnms', 'angle', 'normVx', 'normVy', 'RHO', ...
        'total_ActiveMovie', 'total_AVx', 'total_AVy', '-v7.3');
    
    clear normVx normVy total_ActiveMovie
    
    save([fnm(1:end-4), '_dataSummary.mat'], 'fnms', 'p_Interval1', 'm_Interval1', ...
    'validId', 'durations', 'diameters', 'roiCentr', 'roiArea', 'angle', 'RHO', 'm_p_Duration', 'boundBox', 'pixel', '-v7.3');
end
save([fnm(1:6), '_dataSummary.mat'], 'fnms', 'p_Interval1', 'm_Interval1', ...
    'validId', 'durations', 'diameters', 'roiCentr', 'roiArea', 'angle', 'RHO', 'm_p_Duration', 'boundBox', 'pixel', '-v7.3');


