function [normVx, normVy] = computeFlowField_normalized_xx(imgall, sz, resizeRatio)
% computes optic flow field for imgall, vectors normalized before summing
% up, better for computing the overall directional bias

normVx = [];
normVy = [];

for f = 1:sz(3)-1
    img1 = imgall(:, :, f);
    img2 = imgall(:, :, f+1);
    [Vx, Vy, ~] = opticalFlow(img1, img2);
    Vx = imresize(Vx, resizeRatio, 'bilinear');
    Vy = imresize(Vy, resizeRatio, 'bilinear');

    % Normalize the lengths of the arrows
    mag = sqrt(Vx.^2 + Vy.^2);
    normVx(:, :, f) = Vx ./ mag;
    normVy(:, :, f) = Vy ./ mag;

    id = mag > 0;
    normVx(:, :, f) = normVx(:, :, f) .* id;
    normVy(:, :, f) = normVy(:, :, f) .* id;
end