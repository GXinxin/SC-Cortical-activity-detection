% read imageJ masks and save into .mat format
% clear; clc;

% cd 'Z:\Xinxin\Visual\FRMD7\191126_emxTra2bFrmd7_p10_6.7g_109x'
% cd 'Z:\Xinxin\Visual\JamB\190603_JamBai162_p8_3.9g'
% cd 'Z:\Xinxin\Visual\JamB\190605_2_JamBai162_p10_3.6g'
% cd 'Z:\Xinxin\Visual\JamB\190605_JamBai162_p10_4.2g'
% cd 'Z:\Xinxin\Visual\JamB\190710_JamB_p10_no_3.5g_146x'
% cd 'Z:\Xinxin\Visual\JamB\191016_2_jambAi162_p10_3.2g_124xOriginalMag'
% cd 'Z:\Xinxin\Visual\JamB\191016_jambAi162_p10_3.5g_126xOriginalMag'
% cd 'Z:\Xinxin\Visual\Cdh6\191024_2_cdh_p10_2.2g_126x'
cd 'Z:\Xinxin\Visual\Cdh6\191024_cdh_p10_2.3g_123x'



[filename, pathname] = uigetfile('Choose data file to open');
ROI_original = ReadImageJROI([pathname filename]);


imgSz_y = input('enter the number of rows: ');
imgSz_x = input('enter the number of columns: ');
sz = [imgSz_y, imgSz_x];

fn = input('input .mat mask file name:', 's');


% read ROIs
if length(ROI_original) == 1
    ROI{1} = ROI_original;
else
    ROI = ROI_original;    
end


for i = 1:length(ROI)
    if isfield(ROI{i}, 'mnCoordinates')
        roi{i} = poly2mask(ROI{i}.mnCoordinates(:, 1), ROI{i}.mnCoordinates (:, 2), sz(1), sz(2));   
        roiPolygon{i} = ROI{i}.mnCoordinates;
        roiName{i} = ROI{i}.strName;
    else
        roi{i} = zeros(sz(1), sz(2));
        roi{i}(ROI{i}.vnRectBounds(1):ROI{i}.vnRectBounds(3), ROI{i}.vnRectBounds(2):ROI{i}.vnRectBounds(4)) = 1;
        roiName{i} = ROI{i}.strName;
    end
    roiName{i} = ROI{i}.strName;
end


save([fn, '.mat'], 'roi', 'roiPolygon', 'roiName')
