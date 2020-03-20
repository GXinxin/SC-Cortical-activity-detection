% compute roi event interval for SC and/or V1 recordings
% corrected angle from 032216 version
% add mean active duration for each voxel
% need roi contours, and several movies

% changed interval definition: 
% (previous) on - off
% (new) inter-eventCenter-interval
% 06/22/16


% short version for threshold testing


clear; clc;

% cd 'Z:\Xinxin\Visual_nas6\180202_rxG6_p4_3.5g_female_50x\preprocessed_SC'
% sz = [181 271];
% cd 'Z:\Xinxin\Visual\160505_p3_2.5g_male_52x_relay_switch\preprocessed'
% sz = [218 330];
% cd 'Z:\Xinxin\Visual\160505_p3_2.7g_female_relay\preprocessed'
% sz = [210 330];
% % cd 'Z:\Xinxin\Visual\160508_p4_3.3g_female_49x\preprocessed'
% sz = [191 261];
% cd 'Z:\Xinxin\Visual\160507_p3_female_2.8g_49x\preprocessed'
% sz = [192 292];
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\P10\190616_chatiDTR_p10_5.5g_134x\preprocessed'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\P10\190616_chatiDTR_p10_lowDT_6.3g_63x_upstairs\preprocessed'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\P10\190616_chatiDTR_lowDT_p10_5.5g_116x_oldRig\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\chronicGabazine\190615_emxTra2b_p13_19g_101x_oldRig\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\chronicGabazine\190625_emxTra2b_chronicGabazine_p14_8.2g_102x\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\chronicGabazine\190706_2_EmxTra2b_GabazineChronicInj_p14_63x_upstairs\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\chronicGabazine\190706_EmxTra2b_GabazineChronicInj_p14_63x_upstairs\preprocessed'
% cd 'Z:\Xinxin\Visual_nas6\180322_snap25_p11_6.4g_female_pharm\preprocessed'
% cd 'Z:\Xinxin\Visual_nas6\180325_p9_snap25_5.8g_female_pharm\preprocessed'
% cd 'Z:\Xinxin\Visual_nas6\180402_2_snap25_p11_7.1g_male_upstairs_pharm\preprocessed'
% cd 'Z:\Xinxin\Visual_nas6\180402_snap25_p10_5.9g_female_upstairs_pharm\preprocessed'
% cd 'Z:\Xinxin\Visual_nas6\190711_2_ctrl_chatiDTR_p10_5.6g_upstairs_67.9x_gabazine\preprocessed'
% cd 'Z:\Xinxin\Visual_nas6\190711_ctrl_chatiDTR_p10_5.9g_122x_oldRig\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\chronicGabazine\190712_tra2bSnap25_p10_7.0g_chronicGabazine_108x\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\chronicGabazine\190718_tra2bSnap25_p11_6.9g_chronicGabazine_113x\preprocessed'
% cd 'Z:\Xinxin\Visual\MFA\190718_snap25_p9_10_7.0g_152x\preprocessed'
% cd 'Z:\Xinxin\Visual\MFA\190719_snap25_p11_7.3g_146x\preprocessed'
% cd 'Z:\Xinxin\Visual\MFA\190719_snap25_p11_7.3g_upstairs_57.3x\preprocessed'
% cd 'Z:\Xinxin\Visual\MFA\190720_snap25_p11_5.7g_upstairs_60.6x\preprocessed'
% cd 'Z:\Xinxin\Visual\MFA\190723_p11_snap25_3.2g_upstairs_54.5x\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\190825_tra2bSnap25_p10_6.4g_oldRig_86x\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\190903_tra2bSnap25_p10_8.8g_new71x\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\190904_emxTra2b_p10_5.9g_new71x\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\190910_emxTra2b_p11_6.7g_new71x\preprocessed'
% cd 'Z:\Xinxin\Visual\FRMD7\190928_emxTra2bFrmd7_p10_8.9g_oldRig_82x\preprocessed'
% cd 'E:\Lab\Data\ChAT\Tra2b\170321_emxTra2b_p11_7.0g_male_95x_3rdRig'
% cd 'Z:\Xinxin\Visual\Tra2b\170217_tra2bsnap25_p9_6.0g_male_101x_3rdRig'
% cd 'Z:\Xinxin\Visual\Tra2b\170606_emxTra2b_p9_5.6g_male_gcamp_chrimson_3rdRig\preprocessed'
% cd 'Z:\Xinxin\Visual\Tra2b\170405_emxTra2b_p9_late_6.4g_male_48x_upstairs'
% cd 'Z:\Xinxin\Visual\FRMD7\191126_emxTra2bFrmd7_p10_6.7g_109x\preprocessed'
% cd 'Z:\Xinxin\Visual\FRMD7\191125_emxTra2bFrmd7_p10_5.6g_102x\preprocessed'
% cd 'Z:\Xinxin\Visual\JamB\190603_JamBai162_p8_3.9g\preprocessed'
% cd 'Z:\Xinxin\Visual\JamB\190605_2_JamBai162_p10_3.6g\preprocessed'
% cd 'Z:\Xinxin\Visual\JamB\190605_JamBai162_p10_4.2g\preprocessed'
% cd 'Z:\Xinxin\Visual\JamB\190710_JamB_p10_no_3.5g_146x\preprocessed'
% cd 'Z:\Xinxin\Visual\JamB\191016_2_jambAi162_p10_3.2g_124xOriginalMag\preprocessed'
% cd 'Z:\Xinxin\Visual\JamB\191016_jambAi162_p10_3.5g_126xOriginalMag\preprocessed'
% cd 'Z:\Xinxin\Visual\Cdh6\191024_2_cdh_p10_2.2g_126x\preprocessed'
cd 'Z:\Xinxin\Visual\Cdh6\191024_cdh_p10_2.3g_123x\preprocessed'


% sz = [250 256];
% sz = [270 320];
sz = [256 250];



isSVD = 0;


filelist = readtext('files.txt', ' ');
fnms = filelist(:, 1);
mask_fnms = filelist(:, 2);
if isSVD
    svd_fnms = filelist(:, 3);
end

no_movies = length(fnms);
downSampleRatio = 1;
maskDownSample = 0.5;
sigma = 1; % sigma for gaussian smoothing
% th = [0.75 1 1.5 2 3];
th = [1 1.5 2 3 5]; % segmentation threshold
frameRate = 10; % Hz
% sz = zeros(1, 3);
% sz(1:2) = [270 320];
% sz = [256 250];
% sz = [251 401];
% sz = [512 500];

isTopHat = 1;
diskSz = 20;



for n = [1 7]
    
    
    % load masks for SCs
    load(mask_fnms{n});

    % load movie (preprocessed)
    fnm = fnms{n};
    load(fnm)
    
    if ~exist('ddiff', 'var')
        imgall = dA(:, 100:2000);
        clear dA
    else
        imgall = ddiff(:, 100:2000);
        clear ddiff
    end
    
    imgall = double(imgall);
    imgall(isnan(imgall)) = 0;
    
    
    sz(3) = size(imgall, 2);
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
    imgall = zscore(reshape(imgall, 1, sz(1)*sz(2)*sz(3)));
    imgall = reshape(imgall, sz(1), sz(2), sz(3));
    
    
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




    for t = 1:length(th)
        thresh = th(t);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Segmentation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for r = 1:length(roi)
    %     for r = 2

            savefn = [fnm(1:end-4), '_', num2str(r)];

            if r == length(roi) + 1
                maskMatrix = totalMask;
    %             maskMatrix = imresize(totalMask, maskDownSample, 'bilinear');
            else
                maskMatrix = roi{r};
                if adjust
                    maskMatrix = imresize(roi{r}, .5, 'bilinear');
                end
    %             maskMatrix = imresize(roi{r}, maskDownSample, 'bilinear');
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
            temp_subMov = zscore(temp_subMov');
            temp_subMov = temp_subMov';
            subMov = zeros(size(subMov));
            subMov(maskId, :) = temp_subMov;
    %         if r == length(roi)
    %             thresh = 5;
    %         end
%             activeMov = (subMov > thresh) + (subMov < -thresh);
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
            area = vertcat(STATS.Area);


            validId{r, n} = (dura > 3 & dia > 5); 

            durations{r, n} = dura(validId{r, n});
            diameters{r, n} = dia(validId{r, n});
            roiCentr{r, n} = vertcat(STATS(validId{r, n}).Centroid);
            roiArea{r, n} = [STATS(validId{r, n}).Area];  
            boundBox{r, n} = [STATS(validId{r, n}).BoundingBox];

            valid{r, n} = find(validId{r, n} > 0);



            % reconstruct binary mov based on valid connected components
            valid_mov = subMov;
            valid_activeMov{r} = activeMov;
            badId = find(validId{r, n} == 0);
            for id = 1:length(badId)
                removePixel = CC.PixelIdxList{badId(id)};
                valid_mov(removePixel) = 0;
                valid_activeMov{r}(removePixel) = 0;
            end

            valid_activeMov_down = imresize(valid_activeMov{r}, 0.25, 'bilinear');
            tmpsz = size(valid_activeMov_down);
            valid_activeMov_down = reshape(valid_activeMov_down, tmpsz(1)*tmpsz(2), tmpsz(3));
            
            if r == 1
                total_ActiveMovie = valid_activeMov{r};
            else
                total_ActiveMovie = total_ActiveMovie + valid_activeMov{r};
            end
                          
        end 
        
        
        
        % create segmented movie  
        total_ActiveMovie = total_ActiveMovie > 0;
        for fr = 1:sz(3)
            [I2, map2] = gray2ind(total_ActiveMovie(:,:,fr), 8); %figure; imshow(I2,map)
            F(fr) = im2frame(I2,map2);  %setup the binary segmented mask movie
        end
        fnm3 = [fnm(1:end-4), '_mask_th', num2str(thresh), '.avi'];
        writeMovie_xx(F, fnm3, 0);
        


    end 
end
% save([fnm(1:6), '_dataSummary.mat'], 'fnms', 'p_Interval1', 'm_Interval1', 'p_Interval2', 'm_Interval2', ...
%     'validId', 'durations', 'diameters', 'roiCentr', 'roiArea', 'angle', 'RHO', 'm_p_Duration', 'n_Vx', 'n_Vy');


