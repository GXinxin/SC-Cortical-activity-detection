function [AVx, AVy] = computeFlowField_xx(imgall, sz)
% compute flow field

for fr = 1:sz(3)-1; %option:parfor
    img1 = imgall(:, :, fr);
    img2 = imgall(:, :, fr+1);
    [Vx, Vy, ~] = opticalFlow(img1, img2);
    AVx(:, :, fr) = Vx;
    AVy(:, :, fr) = Vy;
end
AVx(:, :, sz(3)) = AVx(:, :, end);
AVy(:, :, sz(3)) = AVy(:, :, end);