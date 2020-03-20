
clear; clc;



% cd 'Z:\Xinxin\Visual\FRMD7\180702_FRMD7_P9_180704_reimage\results'
% cd 'Z:\Xinxin\Visual\FRMD7\180703_2_P10_FRMD7_retinal_5.1g\results'
% cd 'Z:\Xinxin\Visual\FRMD7\180702_P9_FRMD7_sinus+retinal_4.6g\result'
% cd 'Z:\Xinxin\Visual\FRMD7\180702_2_P9_FRMD7_ctrl_retinal_4.6g\result'
% cd 'Z:\Xinxin\Visual\FRMD7\180703_P10_FRMD7_retinal_5.6g\results'


% cd 'E:\Lab\Data\ChAT\Tra2b\171219_tra2bSnap25_8.4g_male_p9_98x\results'

% cd 'E:\Lab\Data\ChAT\Tra2b\170529_tra2bSnap25_p3_male_2.6g_3rdRig\results'
% cd 'E:\Lab\Data\ChAT\Tra2b\170918_tra2bSnap25_p4_2.3g_male\results'
% cd 'E:\Lab\Data\ChAT\Tra2b\170712_tra2bSnap25_p3_3.6g_female_76x\results'
% colorId = 1;

% cd 'E:\Lab\Data\ChAT\Tra2b\170311_tra2bSnap25_p6_4.0g_female_102x'
% cd 'E:\Lab\Data\ChAT\Tra2b\170403_emxTra2b_p7_4.8g_male_48x_upstairs\results'
% cd 'E:\Lab\Data\ChAT\Tra2b\170215_tra2bsnap25_p7_female_4.4g_48x\results'
% colorId = 2;

% cd 'E:\Lab\Data\ChAT\Tra2b\170130_emxTra2b_p8_6.1g_male_97x_oldRig'
% cd 'E:\Lab\Data\ChAT\Tra2b\170131_emxTra2b_p9_male_6.6g_49.5x_upstair\results'
% cd 'E:\Lab\Data\ChAT\Tra2b\170203_emxTra2b_p9_7g_male_50x_upstairs'
% cd 'E:\Lab\Data\ChAT\Tra2b\170202_emxTra2b_p8_6.8g_male_50x_upstairs\results'
% colorId = 3;

% cd 'E:\Lab\Data\ChAT\Tra2b\170407_emxTra2b_p11_9.1g_female_51x_upstairs\results'
% cd 'E:\Lab\Data\ChAT\Tra2b\170218_emxTra2b_p10_7g_male_53x_sponUpstair_light3rdRig\results'
% cd 'E:\Lab\Data\ChAT\Tra2b\170223_emxTra2b_p10_8.7g_female_125x_3rdRig'
% cd 'E:\Lab\Data\ChAT\Tra2b\170321_emxTra2b_p11_7.0g_male_95x_3rdRig\results'
% colorId = 4;

% cd 'E:\Lab\Data\ChAT\Tra2b\170221_emxTra2b_p13_8.7g_female_115x_3rdRig'
% cd 'E:\Lab\Data\ChAT\Tra2b\170425_tra2bsnap25_8.2g_male_122x\results'
% cd 'E:\Lab\Data\ChAT\Tra2b\170503_emxTra2b_p14_10g_male_118.x_3rdRig\results'
% cd 'E:\Lab\Data\ChAT\Tra2b\170517_emxTra2b_male_p14_8.9g_107x_3rdRig\results'
% colorId = 5;

% cd 'E:\Lab\Data\ChAT\151223\results'

% cd 'E:\Lab\Data\ChAT\Tra2b\170223_emxTra2b_p10_8.7g_female_125x_3rdRig'
% cd 'E:\Lab\Data\ChAT\151223'
% cd 'E:\Lab\Data\ChAT\160203\results'
% colorId = 3;


cmp = colormap(parula(8));
cmp = cmp(2:8, :);
cmapRange = {[0 4 1], [0 3 1], [0 2 1], [0 6 4], [0 5 4]};


% cd 'Z:\Xinxin\Visual\FRMD7\P10_11\181005_2_FRMD7_p10_11_6.3g_female\results'
% cd 'Z:\Xinxin\Visual\FRMD7\P10_11\181005_FRMD7_p10_11_7.2g_138x\results'
% cd 'Z:\Xinxin\Visual\FRMD7\ctrl\P10_11\181009_01_ctrl_FRMD7_p10_6.0g\results'
% cd 'Z:\Xinxin\Visual\FRMD7\ctrl\P10_11\181003_ctrl_frmd7_5.6g_p11_upstairs_50x\results'
% cd 'Z:\Xinxin\Visual\FRMD7\ctrl\P10_11\180916_3_ctrl_FRMD7_p11_5.4g_upstairs_63x\results'
% cd 'Z:\Xinxin\Visual\FRMD7\ctrl\P10_11\180703_3_P10_FRMD7_ctrl_retinal_5.6g\results'
% cd 'Z:\Xinxin\Visual\FRMD7\ctrl\P10_11\180724_6_litter2_ctrl_frmd7_p10_5.2g_132x\results'
% cd 'Z:\Xinxin\Visual\FRMD7\ctrl\P10_11\180723_3_ctrl_frmd7_p10_2.5g_132x\results'
% cd 'Z:\Xinxin\Visual\FRMD7\ctrl\P10_11\180725_7_litter2_ctrl_frmd7_p11_5.2g_152x\results'
% cd 'Z:\Xinxin\Visual\FRMD7\ctrl\P10_11\180916_3_ctrl_FRMD7_p11_5.4g_upstairs_63x\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\P10\180925_chatiDTR_p10_5.4g_140x\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\P10\181018_2_chatiDTR_p10_6.4g_150x\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\P10\181018_chatiDTR_p10_6.4g_146x\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\ctrl\DTinj\180928_inj_DTinj_p10_5.0g_upstair_63x\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\ctrl\DTinj\181218_DTinj_snap25_p10_5.4g\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\ctrl\PBSinj\180928_2_inj_PBSinj_p10_5.2g_upstair_50x\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\ctrl\PBSinj\181028_ctrlPBSinj_chatiDTR_p10_5.2g_upstiars_63x\results'
% cd 'Z:\Xinxin\Visual\chatVgat\190222_chatVgat_p10_7.6g_125x_no1\results'
% cd 'Z:\Xinxin\Visual\chatVgat\190222_2_chatVgat_ctrl_p10_7.1g_125x_no2\results'
% cd 'Z:\Xinxin\Visual\chatVgat\190223_chatVgat_p11_5.3g_120x\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\ctrl\DTinj\180928_inj_DTinj_p10_5.0g_upstair_63x\results'
% cd 'Z:\Xinxin\Visual\ChAT_iDTR\ctrl\DTinj\180927_3_snap25_p9_DTinj_4.8g_142x\results'
% cd 'Z:\Xinxin\Visual\170424_chatsnap25_p10_female_5.2g_65x\results'
% cd 'Z:\Xinxin\Visual\170215_chatsnap25_chrimson_p10_5.2g_female_oldRig\results'
% cd 'Z:\Xinxin\Visual\170216_chatsnap25_chrimson_p11_male_5.7g_3rdRig\results'

fd_list = readtext('E:\Lab\Data\ChAT\Tra2b\AgeGroupAnalysis\new\animalList.txt', ' ');
% fd_list = readtext('E:\Lab\Data\ChAT\iDTR\addAnimalList.txt', ' ');
% fd_list = readtext('E:\Lab\Data\ChAT\Tra2b\pharmacology\animalList.txt', ' ');
% fd_list = readtext('E:\Lab\Data\ChAT\FRMD7\addAnimalList.txt', ' ');
% fd_list = readtext('E:\Lab\Data\ChAT\pharm_nonCortexless\animalList.txt', ' ');


save_path = 'E:\Lab\Data\ChAT\Tra2b\AgeGroupAnalysis\new\OpticFlowPlots\';


for fd = 25:size(fd_list, 1)
    
    clear  normVx_mean  normVy_mean
    
    cd(fd_list{fd, 1})
    colorId = fd_list{fd, 2};
    filelist = dir(fullfile('*1*opticFlow.mat'));
    
% %     if fd == 8
% %         startId = 9;
% %     else
% %         startId = 7;
% %     end
% %     endId = length(filelist);
%     startId = 1;
%     endId = 6;
%     if fd == 11
%         endId = 7;
%     end
% 
%     
%     tag = 1;
    
    startId = fd_list{fd, 3};
    endId = fd_list{fd, 4};
    tag = fd_list{fd, 5};

    % maskFn = dir(fullfile('v1_movie*.mat'));
    maskFn = dir(fullfile('*roiSet*.mat'));
    load(maskFn(tag).name)
    
    % roiId = [1 2];
    totalMask = zeros(size(roi{1}));
    for r = 1:length(roi)
        totalMask = totalMask + roi{r};
    end
    totalMask_large = imresize(totalMask, .5, 'nearest');
    totalMask = imresize(totalMask, .25, 'nearest');
    
    
    % to match colors...
    colorRange = cmapRange{colorId}(1:2);
    totalMask_large = totalMask_large * cmapRange{colorId}(3);
    totalMask = totalMask * cmapRange{colorId}(3);
    
    
    
    index = startId:endId;
    
%     summaryfn = [];
    summaryfn = dir(fullfile('summarized_opticFlow.mat'));
    if ~isempty(summaryfn)
        load(summaryfn(1).name)
    else
        
        % for i = startId : endId
        for i = 1:length(index)
            i
            % for i = 1:min(5, length(filelist))
            
            load(filelist(index(i)).name)
%             fnm = fnms{i};
            
            %     if ~exist('total_AVx')
            %         total_AVx = AVx;
            %         total_AVy = AVy;
            %     end
            %
            
            normVx_mean(:, :, i) = mean(normVx, 3);
            normVy_mean(:, :, i) = mean(normVy, 3);
            %     AVx_mean(:, :, i) = mean(total_AVx, 3);
            %     AVy_mean(:, :, i) = mean(total_AVy, 3);
            
            %     h = figure; imagesc(totalMask > 0); colormap parula(7); caxis(colorRange); hold on
            %     quiver(imresize(normVx_mean(:, :, i), .5, 'bilinear') .* totalMask, imresize(normVy_mean(:, :, i), .5, 'bilinear')...
            %         .* totalMask);
            %     axis image
            %     title(fnm(1:end-4))
            %     set(h, 'Position', [0, 0, 1200, 900]);
            %     h.PaperPositionMode = 'auto';
            %     print([fnm(1:end-4), '_quiver_Norm_S'], '-dpng', '-r0')
            
            
            %     h = figure; imagesc(totalMask > 0); colormap parula(7); caxis(colorRange); hold on
            %     quiver(imresize(AVx_mean(:, :, i), .5, 'bilinear') .* totalMask, imresize(AVy_mean(:, :, i), .5, 'bilinear')...
            %         .* totalMask);
            %     axis image
            %     title(fnm(1:end-4))
            %     set(h, 'Position', [0, 0, 1200, 900]);
            %     h.PaperPositionMode = 'auto';
            %     print([fnm(1:end-4), '_quiver_noNorm_S'], '-dpng', '-r0')
            
            
            
            % compute divergence
            sz = size(normVx_mean);
            
            %     [x, y] = meshgrid(1:sz(2), 1:sz(1));
            %     div(:, :, i) = divergence(squeeze(normVx_mean(:, :, i)), squeeze(normVy_mean(:, :, i)));
            %     d = div(:, :, i) .* totalMask_large;
            %
            %     curls_avg(:, :, i) = curl(squeeze(normVx_mean(:, :, i)), squeeze(normVy_mean(:, :, i)));
            %     c = curls_avg(:, :, i) .* totalMask_large;
            %
            
            %     h = figure;
            %     set(h, 'Position', [0, 0, 1200, 900]);
            %
            %     subplot(1, 3, 1);
            %     quiver(imresize(normVx_mean(:, :, i), .5, 'bilinear') .* totalMask, imresize(normVy_mean(:, :, i), .5, 'bilinear')...
            %     .* totalMask); axis image
            %     set(gca, 'Ydir', 'reverse')
            %
            %     subplot(1, 3, 2); imagesc(d); axis image; caxis([-0.05, 0.05])
            %     subplot(1, 3, 3); imagesc(c); axis image
            %     saveas(h, [fnm(1:end-4), '_sourceCurl.png'])
            
            %     h = figure; imagesc((imresize(totalMask, .5) > 0) * cmapRange{colorId}(3)); colormap parula(7); caxis(colorRange); hold on
            %     quiver(imresize(normVx_mean(:, :, i), .25, 'bilinear') .* imresize(totalMask, .5), imresize(normVy_mean(:, :, i), .25, 'bilinear')...
            %         .* imresize(totalMask, .5), 'Color', 'k', 'MaxHeadSize', 5);
            %     axis image
            %     title('AVG Quiver Norm')
            %     set(h, 'Position', [0, 0, 1200, 900]);
            %     saveas(h, 'AVG_Quiver_Norm_S.png')
            %     h.PaperPositionMode = 'auto';
            %     print('AVG_Quiver_Norm_S', '-dpng', '-r0')
        end
    end
    
    summaryfn = dir(fullfile('*_dataSummary.mat'));
    load(summaryfn(1).name)
    
    for r = 1:length(roi)
        h = figure;
        a = cell2mat(angle(r, index));
        polarhistogram(2*pi - a, 20, 'normalization', 'probability', 'FaceColor', cmp(colorId, :));
        %             set(gca,'YDir','reverse'); % 90 degree is moving downward
        ax = gca;
        ax.LineWidth = 2;      
        title(['  Roi:', num2str(r)]);
        saveas(h, [num2str(tag), '_rosePlot', num2str(r), '.png'])
    end
    
    
    h = figure; imagesc((totalMask > 0) * cmapRange{colorId}(3)); colormap parula(7); caxis(colorRange); hold on
    quiver(imresize(mean(normVx_mean, 3), .5, 'bilinear') .* totalMask, imresize(mean(normVy_mean, 3), .5, 'bilinear')...
        .* totalMask, 'Color', 'k', 'MaxHeadSize', 5);
    axis image
    title('AVG Quiver Norm')
    set(h, 'Position', [0, 0, 1200, 900]);
    saveas(h, [num2str(tag), 'AVG_Quiver_Norm_S.png'])
%     saveas(h, [save_path, num2str(fd), 'AVG_Quiver_Norm_S.png'])
    
    
    totalMask0 = imresize(totalMask, .5, 'bilinear');
    h = figure; set(h, 'position', [200 100 700 600])
    imagesc((totalMask0 > 0) * cmapRange{colorId}(3)); colormap parula(7); caxis(colorRange); hold on
    hq = quiver(imresize(mean(normVx_mean, 3), .25, 'bilinear') .* totalMask0, imresize(mean(normVy_mean, 3), .25, 'bilinear')...
        .* totalMask0, 'Color', 'k', 'MaxHeadSize', 0.5, 'LineWidth', 0.1, 'AutoScale','on', 'AutoScaleFactor', 1);
    axis image
    U = hq.UData;
    V = hq.VData;
    X = hq.XData;
    Y = hq.YData;
    U = U';
    V = V';
    
    h = figure; set(h, 'position', [200 100 700 600])
    imagesc((totalMask0 > 0) * cmapRange{colorId}(3)); colormap parula(7); caxis(colorRange);
    set(gca,'YDir','normal');hold on
    headLength = 7;
    LineLength = 2;
    vector_length = sqrt(U.^2 + V.^2);
    max_length = max(vector_length(:));
    for ii = 1:length(X)
        for ij = 1:length(Y)
            head_L = headLength * (0.7 * sqrt(U(ii, ij)^2 + V(ii, ij)^2) / max_length + 0.3);
            head_W = head_L;
            if abs(U(ii,ij) * V(ii,ij)) > 0
                ah = annotation('arrow',...
                    'headStyle','cback1','HeadLength',head_L,'HeadWidth',head_W);
                set(ah,'parent',gca);
                set(ah,'position',[X(ii) Y(ij) LineLength*U(ii,ij)/max_length LineLength*V(ii,ij)/max_length]);
                set(ah, 'LineWidth', 1)
            end
        end
    end
    xlim([0 32]); ylim([0 32])
    axis image
    title('AVG Quiver Norm')
    saveas(h, [num2str(tag), 'AVG_Quiver_Norm_SS.png'])
%     saveas(h, [num2str(tag), 'AVG_Quiver_Norm_SS.png'])
    
    % normalized by the same maximum vector length across animals
    h = figure; set(h, 'position', [200 100 700 600])
    imagesc((totalMask0 > 0) * cmapRange{colorId}(3)); colormap parula(7); caxis(colorRange);
    set(gca,'YDir','normal');hold on
    headLength = 7;
    LineLength = 2;
    for ii = 1:length(X)
        for ij = 1:length(Y)
            head_L = headLength * (0.7 * sqrt(U(ii, ij)^2 + V(ii, ij)^2) / 1.5 + 0.3);
            head_W = head_L;
            if abs(U(ii,ij) * V(ii,ij)) > 0
                ah = annotation('arrow',...
                    'headStyle','cback1','HeadLength',head_L,'HeadWidth',head_W);
                set(ah,'parent',gca);
                set(ah,'position',[X(ii) Y(ij) LineLength*U(ii,ij)/1.5 LineLength*V(ii,ij)/1.5]);
                set(ah, 'LineWidth', 1)
            end
        end
    end
    xlim([0 32]); ylim([0 32])
    axis image
    title('AVG Quiver Norm')
%     saveas(h, [num2str(tag), 'AVG_Quiver_Norm_SS_normLegnth.png'])
    
    
    
    
    % downsampled more
    totalMask_ss = imresize(totalMask0, .5, 'bilinear');
    h = figure; set(h, 'position', [200 100 700 600])
    imagesc((totalMask_ss > 0) * cmapRange{colorId}(3)); colormap parula(7); caxis(colorRange); hold on
    hq = quiver(imresize(mean(normVx_mean, 3), .125, 'bilinear') .* totalMask_ss, imresize(mean(normVy_mean, 3), .125, 'bilinear')...
        .* totalMask_ss, 'Color', 'k', 'MaxHeadSize', 0.5, 'LineWidth', 0.1, 'AutoScale','on', 'AutoScaleFactor', 1);
    axis image
    U = hq.UData;
    V = hq.VData;
    X = hq.XData;
    Y = hq.YData;
    U = U';
    V = V';
    
    h = figure;set(h, 'position', [200 100 700 600])
    imagesc((totalMask_ss > 0) * cmapRange{colorId}(3)); colormap parula(7); caxis(colorRange);
    set(gca,'YDir','normal');hold on
    headLength = 10;
    LineLength = 1.5;
    vector_length = sqrt(U.^2 + V.^2);
    max_length = max(vector_length(:));
    for ii = 1:length(X)
        for ij = 1:length(Y)
            head_L = headLength * (0.7 * sqrt(U(ii, ij)^2 + V(ii, ij)^2) / max_length + 0.3);
            head_W = head_L;
            if abs(U(ii,ij) * V(ii,ij)) > 0
                ah = annotation('arrow',...
                    'headStyle','cback1','HeadLength',head_L,'HeadWidth',head_W);
                set(ah,'parent',gca);
                set(ah,'position',[X(ii) Y(ij) LineLength*U(ii,ij)/max_length LineLength*V(ii,ij)/max_length]);
                set(ah, 'LineWidth', 1)
            end
        end
    end
    xlim([0 16]); ylim([0 16])
    axis image
    title('AVG Quiver Norm')
%     saveas(h, [num2str(tag), 'AVG_Quiver_Norm_SSS.png'])
    
    
    % normalized by the same maximum vector length
    h = figure;set(h, 'position', [200 100 700 600])
    imagesc((totalMask_ss > 0) * cmapRange{colorId}(3)); colormap parula(7); caxis(colorRange);
    set(gca,'YDir','normal');hold on
    headLength = 10;
    LineLength = 1.5;
    vector_length = sqrt(U.^2 + V.^2);
    for ii = 1:length(X)
        for ij = 1:length(Y)
            head_L = headLength * (0.7 * sqrt(U(ii, ij)^2 + V(ii, ij)^2) / 1.5 + 0.3);
            head_W = head_L;
            if abs(U(ii,ij) * V(ii,ij)) > 0
                ah = annotation('arrow',...
                    'headStyle','cback1','HeadLength',head_L,'HeadWidth',head_W);
                set(ah,'parent',gca);
                set(ah,'position',[X(ii) Y(ij) LineLength*U(ii,ij)/1.5 LineLength*V(ii,ij)/1.5]);
                set(ah, 'LineWidth', 1)
            end
        end
    end
    xlim([0 16]); ylim([0 16])
    axis image
    title('AVG Quiver Norm')
%     saveas(h, [num2str(tag), 'AVG_Quiver_Norm_SSS_normLegnth.png'])
    
    
    % h = figure;
    % imagesc(mean(div, 3) .* totalMask_large);
    % axis image; caxis([-0.05, 0.05]); colormap jet
    % title('Source/Sink')
    % set(h, 'Position', [0, 0, 1200, 900]);
    % saveas(h, 'Divergence.png')
    % h.PaperPositionMode = 'auto';
    % print('AVG_Quiver_noNorm_S', '-dpng', '-r0')
    
    
    % h = figure;
    % imagesc(mean(curls_avg, 3) .* totalMask_large);
    % axis image; caxis([-0.05, 0.05]); colormap jet
    % title('Curl')
    % set(h, 'Position', [0, 0, 1200, 900]);
    % saveas(h, 'Curl.png')
    
    
    save([save_path, 'animal', num2str(fd), '_summarized_opticFlow.mat'], 'totalMask', 'totalMask_large', 'normVx_mean', 'normVy_mean', 'roi')
    close all
end