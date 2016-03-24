clear; clc;

cd 'D:\Lab\Data\ChAT\In vitro\160226'

downSampleRatio = 0.25;
sigma = 3;
hasROI = 1;
roifn = '160226_4x_flatField.mat';

mov = VideoReader('16022601_4xnoCover.avi');
nFrames=mov.NumberOfFrames;
for i = 1 : nFrames
  a = read(mov, i);
  img1(:, :, i) = rgb2gray(a);
end
% img1 = img1(:, :, 140:end);

mov = VideoReader('16022602_4xCover.avi');
nFrames=mov.NumberOfFrames;
for i = 1 : nFrames
  a = read(mov, i);
  img2(:, :, i) = rgb2gray(a);
end
% 
% mov = VideoReader('MOVIE5_gray.avi');
% nFrames=mov.NumberOfFrames;
% for i = 1 : nFrames
%   a = read(mov, i);
%   img3(:, :, i) = rgb2gray(a);
% end
% img3 = img3(:, :, 1:90);

% imgall = double(cat(3, img1, img2, img3));
% imgall = double(cat(3, img1, img2));
imgall = imgall(1:127, :, :);
sz = size(imgall);

% dF/F
imgall = reshape(imgall, sz(1)*sz(2), sz(3));
avgImg = mean(imgall, 2);
imgall = imgall ./ repmat(avgImg, 1, sz(3)) - 1;
imgall = reshape(imgall, sz(1), sz(2), sz(3));


imgall = imgall(10:end, :, :); % for 160226_4x_noCover
sz = size(imgall);
savefn = '160226_4x_noCover';
% write dF/F
Iarr=mat2gray(imgall);   %scale the whole array
[I2arr, map] = gray2ind(Iarr, 256); %convert the whole array to 8bit indexed		
for fr = 1:sz(3)
    M(fr) = im2frame(I2arr(:,:,fr), map);  %setup the indexed raw dFoF movie
end
fnm0 = [savefn, '_dFF.avi'];
writeMovie(M, fnm0);


% remove slow trend using rolling hat
hat = 30;
se = strel('line', hat, 0);
imgall = reshape(imgall, sz(1)*sz(2), sz(3));
parfor p = 1:sz(1)*sz(2) 
    imgall(p, :) = imtophat(imgall(p, :), se);
end
imgall = reshape(imgall, sz);


% regress out global fluctuation from reference pixels
% ref = squeeze(imgall(10, 10, :)); % for 160225_10x
% ref = squeeze(imgall(1, 66, :)); % for 160225_4x
ref = squeeze(sum(sum(imgall, 1), 2)); % for 160226_4x_noCover
imgall = reshape(imgall, sz(1)*sz(2), sz(3));
parfor i = 1 : sz(1)*sz(2)
    beta(i) = regress(imgall(i, :)', ref);
end
diff = imgall - repmat(beta', 1, sz(3)) .* repmat(ref', sz(1)*sz(2), 1);

parfor i = 1 : sz(1)*sz(2)
    diff2(i, :) = detrend(diff(i, :));
end
diff2 = reshape(diff2, sz(1), sz(2), sz(3));
% figure; plot(squeeze(diff2(50, 50, :)))


% apply gaussian smooth
diff2 = reshape(diff2, sz(1), sz(2), sz(3));
parfor fr = 1:sz(3); 
    diff2(:, :, fr) = gaussSmooth(diff2(:,:,fr), sigma, 'same');
end	




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute flow field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff2 = reshape(diff2, sz(1), sz(2), sz(3));   
[normVx, normVy] = computeFlowField_normalized_xx(diff2, sz, downSampleRatio);

h = figure; 
quiver(mean(normVx, 3), mean(normVy, 3));
title(savefn)
saveas(h, [savefn, '_quiver.png']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

if hasROI
    load(roifn);
    mask = repmat(roi{1}(10:127, :, :), 1, 1, sz(3));
    diff2 = diff2 .* mask;
end

% zscore
thresh = 2;
diff_zscore = zscore(reshape(diff2, 1, sz(1)*sz(2)*sz(3)));
activeMov = diff_zscore > thresh;
diff_zscore = reshape(diff_zscore, sz(1), sz(2), sz(3));
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

validId = (dura > 2 & dia > 10);
durations = dura(validId);
diameters = dia(validId);
roiCentr = vertcat(STATS(validId).Centroid);
roiArea = [STATS(validId).Area];  


% Compute wave direction
[AVx, AVy] = computeFlowField_xx(diff2, sz); 

nDomains = sum(validId);
Vsum(1:nDomains) = struct('theta',[],'rho',[]);
validDomains = CC.PixelIdxList(validId);

for k = 1:nDomains
    [theta, rho]= cart2pol(sum(AVx(validDomains{k})), sum(AVy(validDomains{k})));
    Vsum(k).theta = theta;
    Vsum(k).rho = rho;
end

angle = vertcat(Vsum.theta);
h = figure; rose(angle);
set(gca,'XDir','reverse', 'YDir','reverse'); % 90 degree is moving downwards
title(savefn);
saveas(h, [savefn, '_rosePlot.png'])



% reconstruct binary mov based on valid connected components
valid_activeMov = activeMov;
badId = find(validId == 0);
for id = 1:length(badId)
    removePixel = CC.PixelIdxList{badId(id)};
    valid_activeMov(removePixel) = 0;
end

valid_activeMov_down = imresize(valid_activeMov, downSampleRatio, 'bilinear');
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

%     goodId = pixelDuration{p} > 1;
%     activeOn{p} = activeOn{p}(goodId);
%     activeOff{p} = activeOff{p}(goodId);
%     pixelDuration{p} = pixelDuration{p}(goodId);

    pixelInterval{p} = activeOn{p}(2:end) - activeOff{p}(1:end-1);
    meanInterval(p) = mean(pixelInterval{p});
end

meanInterval(isnan(meanInterval)) = 0;
meanInterval(meanInterval > 600) = 0; % remove super long intervals for a better range of colorbar

meanInterval = reshape(meanInterval, tmpsz(1), tmpsz(2));
h = figure; imagesc(meanInterval); colorbar; colormap jet
title(savefn);
saveas(h, [savefn, '_interval.png']);
