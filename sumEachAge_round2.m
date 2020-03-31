% Integrate the summary data from *dataSummary files, modified to count
% each hemisphere separately, 10/30/18

clear; clc;
 
path = 'E:\Lab\Data\ChAT\Tra2b\AgeGroupAnalysis\new';
ageGroup = {'p3_4', 'p6_7', 'p8_9', 'p10_11', 'p13_14', 'p16_17', 'p28_29'};
cmp = colormap(parula(8));
cmp = cmp(2:8, :);

frameRate = 10;

excludeId = {zeros(3, 2), zeros(3, 2), [0 0; 0 0; 1 0], zeros(4, 2), zeros(4, 2), [1 0; 0 1; 0 0], [1 0; 0 1; 0 0; 0 0]};



for ag = 1 : 7
    
    individual_var{ag} = [];
    individual_circ_mean{ag} = [];
    badId = excludeId{ag};
    
    clear pixel intervals1 total_p_duration totalAngle total_diameter total_duration...
        total_area regionFreq interval_median1 interval_mean1 interval_std1 totalInterval1...
        center_dist total_center_dist total_center_on total_center_stepDist_sum total_center_speed...
        total_center_stepSpeed total_max_area center_on center_stepDist_sum center_speed...
        center_stepSpeed max_area
    
    data_path = [path, '\', ageGroup{ag}];
    cd(data_path)
    filelist = dir(fullfile('*dataSummary*'));
    
    rig3 = 40 * 23.36 * 2; % converted from measurements of the cameras, multiplied by 2 for processed downsampled data
    rig_up = 20 * 20.4 * 2;
    rig_old = 40 * 16.92 * 2;
    
    % magnification for all animals
    mag = {[rig3 / 79, rig3 / 80, rig3 / 88], [rig_up / 48, rig3 / 102, rig_up / 48], ...
        [rig_old / 97, rig_old / 112, rig_up / 50, rig_up / 50], [rig_up / 53, rig3 / 125, rig3 / 95, rig_up / 51], ...
        [rig3 / 115, rig3 / 122, rig3 / 118, rig3 / 107], ...
        [rig3 / 122, rig3 / 112, rig3 / 105], ...
        [rig3 / 122, rig3 / 112, rig3/118, rig3 / 107]};
%     mag = {[rig3 / 79, rig3 / 76, rig3 / 88], [rig_up / 48, rig3 / 102, rig_up / 48], ...
%         [rig_old / 97, rig_old / 112, rig_up / 50, rig_up / 50], [rig_up / 53, rig3 / 125, rig3 / 95, rig_up / 51], ...
%         [rig3 / 115, rig3 / 122, rig3 / 118, rig3 / 107], ...
%         [rig3 / 122, rig3 / 112, rig3 / 105], ...
%         [rig3 / 122, rig3 / 112, rig3 / 107]};

    for f = 1:length(filelist)
        fn = filelist(f).name;
        load(fn)
        allfn = fn(1:9); 
        img_sz = 2 * size(m_p_Duration{1, 1});
        
        
        sz = size(diameters);
        if sz(2) > 6
            sz(2) = 6;
        end
        
        addTmp = cell(1, 6);
        if sz(1) == 1 && badId(f, 1) == 1
            p_Interval1 = [addTmp; p_Interval1];
            m_Interval1 = [addTmp; m_Interval1];
            validId = [addTmp; validId];
            durations = [addTmp; durations];
            diameters = [addTmp; diameters];
            roiCentr = [addTmp; roiCentr];
            roiArea = [addTmp; roiArea];
            angle = [addTmp; angle];
            RHO = [addTmp; RHO];
            m_p_Duration = [addTmp; m_p_Duration];
            boundBox = [addTmp; boundBox];
            pixel = [addTmp; pixel];
            
            sz = size(diameters);
        end

        
        for r = 1 : sz(1) % select which ROIs to look at
            
            total_diameter{f}{r} = [];
            total_duration{f}{r} = [];
            total_area{f}{r} = [];
            total_center_dist{f}{r} = [];
            total_center_stepDist_sum{f}{r} = [];
            total_center_speed{f}{r} = [];
            total_center_stepSpeed{f}{r} = [];
            total_max_area{f}{r} = [];
            total_center_on{f}{r} = [];
            

            if badId(f, r) == 0

                for j = 1:sz(2)
                    % overall stats
                    regionFreq{f}(r, j) = length(diameters{r, j}); % get all wave frequencies
                    total_diameter{f}{r} = [total_diameter{f}{r}; diameters{r, j} * mag{ag}(f)];
                    total_duration{f}{r} = [total_duration{f}{r}; durations{r, j} /  frameRate];
%                     total_area{f}{r} = [total_area{f}{r}; roiArea{r, j}' * (mag{ag}(f)^2)]; 
                    total_area{f}{r} = [total_area{f}{r}; roiArea{r, j}' * (mag{ag}(f)^2) ./ durations{r, j}];                


                    if exist('pixel', 'var')
                        for w = 1:size(pixel{r, j}, 2)

                            minT = min(pixel{r, j}{w}(:, 3));
                            maxT = max(pixel{r, j}{w}(:, 3));
                            timeSteps = minT : maxT;
                            for t = 1 : length(timeSteps)
                                p_id{t} = find(pixel{r, j}{w}(:, 3) == timeSteps(t));
                                length_p{f}{r, j}(w, t) = length(p_id{t});
    %                             p_loc = sub2ind(img_sz, pixel{f}{r, j}{w}(p_id{t}, 1), pixel{f}{r, j}{w}(p_id{t}, 2));
    %                             time_img = zeros(img_sz);
    %                             time_img(p_loc) = 1;
                                center_time{f}{r, j}{w}(t, :) = [mean(pixel{r, j}{w}(p_id{t}, 1)), mean(pixel{r, j}{w}(p_id{t}, 2))];
                                if t > 1
                                    center_stepDist{f}{r, j}{w}(t-1) = pdist([center_time{f}{r, j}{w}(t-1, :); center_time{f}{r, j}{w}(t, :)],'euclidean')...
                                        * mag{ag}(f);
                                end
                            end


                            center_on{f}{r, j}(w, :) = center_time{f}{r, j}{w}(1, :);

                            % centroid travel distance
                            center_dist{f}{r, j}(w) = pdist([center_time{f}{r, j}{w}(1, :); center_time{f}{r, j}{w}(end, :)],'euclidean') * mag{ag}(f); 
                            center_stepDist_sum{f}{r, j}(w) = sum(center_stepDist{f}{r, j}{w});

                            % centroid moving speed
                            center_speed{f}{r, j}(w) = center_dist{f}{r, j}(w) / length(timeSteps) * frameRate; % unit: um/s
                            center_stepSpeed{f}{r, j}(w) = mean(center_stepDist{f}{r, j}{w}) * frameRate; % averaged from each time point, unit: um/s

                            % maximum area during each wave
                            max_area{f}{r, j}(w) = max(length_p{f}{r, j}(w, :)) * (mag{ag}(f)^2);

                        end  
                    end
                    
                    total_center_dist{f}{r} = [total_center_dist{f}{r}, center_dist{f}{r, j}];
                    total_center_stepDist_sum{f}{r} = [total_center_stepDist_sum{f}{r}, center_stepDist_sum{f}{r, j}];
                    total_center_speed{f}{r} = [total_center_speed{f}{r}, center_speed{f}{r, j}];
                    total_center_stepSpeed{f}{r} = [total_center_stepSpeed{f}{r}, center_stepSpeed{f}{r, j}];
                    total_max_area{f}{r} = [total_max_area{f}{r}, max_area{f}{r, j}];
                    total_center_on{f}{r} = [total_center_on{f}{r}; center_on{f}{r, j}]; 
                        
                end
            
            
                
                h = figure; 
                set(h, 'visible', 'off')
                scatter(total_diameter{f}{r}*2, total_duration{f}{r}/frameRate)
                xlabel('diameter (pixel)')
                ylabel('duration (s)')
                saveas(h, [allfn, '_size.png'])
            end
        end

        
        
        
        % pixel wave interval
        totalInterval1{f}{r} = zeros(1, img_sz(1)*img_sz(2)/4);
        for p = 1:length(p_Interval1{1, 1})
            for r = 1:size(angle, 1)   
                if badId(f, r) == 0
                    interval = [];
                    for id = 1:sz(2)
                        interval = [interval, p_Interval1{r, id}{p}];
                    end
                    totalInterval1{f}{r}(p) = mean(interval);    
                end
            end
        end


        % pixel active duration
        for r = 1:size(angle, 1)       
            if badId(f, r) == 0
                d = zeros(size(m_p_Duration{r, 1}));
                for id = 1:sz(2)
                    if ~isempty(m_p_Duration{r, id})
                        d = d + m_p_Duration{r, id};
                    end
                end
                total_p_duration{f}{r} = d / sz(2);  
            end
        end



        % wave angle
        for r = 1:size(angle, 1)
            if badId(f, r) == 0
                totalAngle{f}{r} = [];
                for id = 1:sz(2)
                    totalAngle{f}{r} = [totalAngle{f}{r}, angle{r, id}];
                end
                
                savefn = [allfn, '_Roi', num2str(r)];
                
                
                tmp = wrapTo2Pi(angle{r, id}');    
                individual_var{ag} = [individual_var{ag}, circ_var(tmp)];
                individual_circ_mean{ag} = [individual_circ_mean{ag}, circ_mean(tmp)];
                
                
                
                
                h = figure;
                polarhistogram(2*pi - totalAngle{f}{r}, 20, 'normalization', 'probability', 'LineWidth', 1);
                thetaticks(0:45:315);
                rticks([0.05, 0.1, 0.15])
                rticklabels({})
                rlim([0 0.15])
                ax = gca;
                ax.LineWidth = 1;
                title(['Roi:', num2str(r)]);
                saveas(h, [savefn, '_rosePlot.png'])
                
                
                h = figure; set(h, 'visible', 'off')
                imagesc(reshape(totalInterval1{f}{r}/frameRate, img_sz(1)/2, img_sz(2)/2));               
                colorbar; caxis([0 120]); colormap jet
                title(['Roi:', num2str(r)])
                saveas(h, [savefn, '_meanInterval1.png'])
                
                
                h = figure; set(h, 'visible', 'off')
                imagesc(total_p_duration{f}{r});
                colorbar; caxis([0, 0.2]); colormap jet
                title(['Roi:', num2str(r)])
                saveas(h, [savefn, '_duration.png'])
                
                
%                 validP = (totalInterval1{f}{r} > 10) & (totalInterval1{f}{r} < 2000);
%                 intervals1{f}{r, j} = totalInterval1{f}{r}(validP)/frameRate;
%                 interval_mean1(r) = mean(intervals1{f}{r, j});
%                 interval_median1(r) = median(intervals1{f}{r, j});
%                 interval_std1(r) = std(intervals1{f}{r, j})/sqrt(length(intervals1{f}{r, j}));
%                 
%                 %
%                 p_durations{f}{r, j} = total_p_duration{f}{r, j}(validP);
%                 p_durations_mean(r, j) = mean(p_durations{f}{r, j});
%                 %         interval_std(r, j) = std(intervals{f}{r, j});
%                 p_durations_std(r, j) = std(p_durations{f}{r, j})/sqrt(length(p_durations{f}{r, j}));
%                 
                
%                 % interval
%                 h = figure;
%                 set(h, 'visible', 'off')
%                 errorbar(interval_mean1(:, r), interval_std1(:, r));
%                 title(['Roi:', num2str(r)])
%                 %             ylim([30 90])
%                 xlim([0 5])
%                 savefn2 = [allfn, '_Roi', num2str(r)];
%                 saveas(h, [savefn2, '_intervalAVG1.png'])
%                 
%                 
%                 
%                 figure; hist(intervals1{r, j}, 0:15:120);
%                 h = findobj(gca,'Type','patch');
%                 set(h,'EdgeColor','w','facealpha',0.6)
%                 saveas(h, [savefn2, '_intervalHist1.png'])
                
            end
        end
    end
    
    
    save([path, '\', ageGroup{ag}, '_summaryStat_test_010820.mat'], ... 
        'totalInterval1', 'total_p_duration', 'totalAngle', 'total_diameter', 'total_duration', 'total_area', 'regionFreq', ...
        'total_center_dist', 'total_center_on', 'total_center_stepDist_sum', 'total_center_speed', 'total_center_stepSpeed', 'total_max_area', ...
        'center_dist', 'center_on', 'center_stepDist_sum', 'center_speed', 'center_stepSpeed', 'max_area', 'individual_var', 'individual_circ_mean')
    close all
 
end



