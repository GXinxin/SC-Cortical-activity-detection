clear; clc

cd 'E:\Lab\Data\ChAT\Tra2b\AgeGroupAnalysis\new'
% %
%
% ageGroup = {'p3_4', 'p6_7', 'p8_9', 'p10_11', 'p13_14', 'p16_17', 'p28_29'};
ageGroup = {'p3_4', 'p6_7', 'p8_9', 'p10_11', 'p13_14', 'p28_29'};
cmp = colormap(parula(8));
cmp = cmp(2:8, :);
%
% %
% badId = cell(1, 7);
% badId{3} = [2, 1; 4, 1];

% cd 'E:\Lab\Data\ChAT\FRMD7\exp'
% cd 'E:\Lab\Data\ChAT\FRMD7\ctrl'
% ageGroup = {'p9', 'p10'};
% cmp = colormap(gray(6));
% cmp = cmp([4, 2], :);

% cd 'E:\Lab\Data\ChAT\iDTR'
% ageGroup = {'ctrl_noinj', 'ctrl_pbs', 'exp'};
% cmp = colormap(pink(6));
% cmp = cmp([1, 2, 4], :);

% cd 'E:\Lab\Data\ChAT\B2\baylorb2\summary'
% ageGroup = {'p9'};
% cmp = colormap(gray(6));
% cmp = cmp([3], :);

%
% cd 'E:\Lab\Data\ChAT\RxG6'
% ageGroup = {'p3_4', 'p6_7', 'p8_9', 'p10_11'};
% cmp = colormap(parula(8));
% cmp = cmp(2:8, :);


filelist = dir(fullfile('*summaryStat_test_010820.mat'));
fs = 10;

% badId = {zeros(3, 2), [0 0; 0 0; 0 1], [0 0; 0 0; 1 0], zeros(4, 2), zeros(4, 2), [1 0; 0 1; 0 0], [1 0; 0 1; 0 0; 0 0]};
badId = {zeros(3, 2), [0 0; 0 0; 0 1], [0 0; 0 0; 1 0], zeros(4, 2), zeros(4, 2), [1 0; 0 1; 0 0; 0 0]};




for age = 1:length(filelist)
    
    load(filelist(age).name)
    savefn = ['age_', ageGroup{age}];
    
    sz(1) = size(totalAngle, 2);
    sz(2) = 2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANGLE
    % generate combined rose plot for each hemisphere
    for r = 1:sz(2)
        angles{age, r} = [];
        for n = 1:sz(1)
            if badId{age}(n, r) == 0
                angles{age, r} = [angles{age, r}, totalAngle{n}{r}];
            end
        end
        %         h = figure;
        %         polarhistogram(2*pi - angles{age, r}, 20, 'normalization', 'probability');
        %         title(['Roi:', num2str(r)]);
        %         saveas(h, [savefn, 'roi', num2str(r), '_rosePlot.png'])
        %         saveas(h, [savefn, 'roi', num2str(r), '_rosePlot.eps'])
    end
    
    % generate combined rose plot for both hemispheres together, on left SC reference
    both_angles{age} = [angles{age, 1}, 3*pi - angles{age, 2}];
    h = figure;
    polarOutput(age) = polarhistogram(2*pi - both_angles{age}, 20, 'normalization', 'probability');
    %     polarhistogram(2*pi - both_angles{age}, 20, 'normalization', 'probability', 'EdgeColor', cmp(age, :), 'FaceColor', [1 1 1]);
    polarhistogram(2*pi - both_angles{age}, 20, 'normalization', 'probability', 'FaceColor', cmp(age, :));
%     thetaticks(0:45:315);
%     rticks([0.05, 0.1, 0.15])
%     rticklabels({})
%     rlim([0 0.15])
    ax = gca;
    ax.LineWidth = 2;
    title('both');
    set(gca,'FontSize', 16);
%     saveas(h, [savefn, '_both_rosePlot.png'])
%     saveas(h, [savefn, '_both_rosePlot.eps'])
%     
    
    tmp = wrapTo2Pi(both_angles{age}');    
    total_var(age) = circ_var(tmp);
    total_circ_mean(age) = circ_mean(tmp);
    total_circ_std(age) = circ_std(tmp);
    
    
    
    %%%%%%%%%%%%%%%%%%%%
    % AREA
    
    for r = 1:sz(2)
        Areas{age, r} = [];
        Durations{age, r} = [];
        Intervals{age, r} = [];
        for n = 1:sz(1)
            if badId{age}(n, r) == 0
                total_area{n}{r} = total_area{n}{r} ./ total_duration{n}{r};
                medianArea(age, n, r) = median(total_area{n}{r});
                medianArea_max(age, n, r) = median(total_max_area{n}{r});
                Areas{age, r} = [Areas{age, r}; total_area{n}{r}];
                
                medianDuration(age, n, r) = median(total_duration{n}{r});
                Durations{age, r} = [Durations{age, r}; total_duration{n}{r}];
                
                validP = (totalInterval1{n}{r} > 10) & (totalInterval1{n}{r} < 2000);
                intervals1{n, r} = totalInterval1{n}{r}(validP) / fs;
                medianInterval(age, n, r) = median(intervals1{n, r});
                Intervals{age, r} = [Intervals{age, r}, intervals1{n, r}];
                
                median_center_dist(age, n, r) = median(total_center_dist{n}{r});
                median_center_stepDist_sum(age, n, r) = median(total_center_stepDist_sum{n}{r});
                median_center_speed(age, n, r) = median(total_center_speed{n}{r});
                median_center_stepSpeed(age, n, r) = median(total_center_stepSpeed{n}{r});
            end
        end
    end
end


figure;
anova_tmp = []; group_tmp = [];
for i = 1:7%size(ageGroup, 2)
    if i ~=6
        anova_tmp = [anova_tmp; individual_var{i}'];
        group_tmp = [group_tmp; i * ones(length(individual_var{i}), 1)];
    end
end
[CV_anova_p,CV_t,CV_stats] = anova1(anova_tmp, group_tmp);
[CV_c, m] = multcompare(CV_stats);





clear tmp

h = figure;
set(h, 'Position', [0 0 600, 500])
for i = 1:size(ageGroup, 2)
    tmp{i} = squeeze(medianArea(i, :, :));
    tmp{i} = tmp{i}(tmp{i} > 0);
    scatter(i*ones(size(tmp{i})),tmp{i}, [], cmp(i, :), 'LineWidth', 2); hold on
end
set(gca, 'YScale', 'log')

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
boxplot2(tmp, 'Width', .6, 'Colors', cmp, 'Labels', ageGroup)
title('Area')
ylabel('um^2')
set(findobj(gca,'type','line'),'linew',3)
set(gca,'linew',3, 'FontSize', 22)
% ylim([10000, 1000000])
saveas(h, 'groupedArea_test.png')
saveas(h, 'groupedArea_test.eps')






h = figure;
set(h, 'Position', [0 0 600, 500])
anova_tmp = [];
group_tmp = [];
for i = 1:size(ageGroup, 2)
    tmp{i} = squeeze(medianArea_max(i, :, :));
    tmp{i} = tmp{i}(tmp{i} > 0);
    scatter(i*ones(size(tmp{i})),tmp{i}, [], cmp(i, :), 'LineWidth', 2); hold on
    anova_tmp = [anova_tmp; tmp{i}];
    group_tmp = [group_tmp; i * ones(length(tmp{i}), 1)];
end
set(gca, 'YScale', 'log')

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
boxplot2(tmp, 'Width', .6, 'Colors', cmp, 'Labels', ageGroup)
box off
title('Area Peak')
ylabel('um^2')
set(findobj(gca,'type','line'),'linew',3)
set(gca,'linew',3, 'FontSize', 22)
ylim([5000, 320000])
set(gca,'TickLength',[0.015, 0.001])
saveas(h, 'groupedArea_peak_test.png')
saveas(h, 'groupedArea_peak_test.eps')

% comb = nchoosek(1:6, 2);
% star_area = zeros(1, size(comb, 1));
% for c = 1 : size(comb, 1)
%     [h_area(c), p_area(c)] = ttest2(tmp{comb(c, 1)}, tmp{comb(c, 2)});
% end
% star_area(p_area < 0.05) = 1;
% star_area(p_area < 0.01) = 2;
% star_area(p_area < 0.001) = 3;

[area_anova_p,area_t,area_stats] = anova1(anova_tmp, group_tmp, 'off');
[area_c, m] = multcompare(area_stats);
area_t
for i = 1 : length(tmp)
    aa(1, i) = mean(tmp{i});
    aa(2, i) = std(tmp{i})/sqrt(length(tmp{i}));
end



h = figure;
set(h, 'Position', [0 0 600 500])
for i = 1:size(ageGroup, 2)
    tmp{i} = squeeze(median_center_dist(i, :, :));
    tmp{i} = tmp{i}(tmp{i} > 0);
    scatter(i*ones(size(tmp{i})),tmp{i}, [], cmp(i, :), 'LineWidth', 2); hold on
end
set(gca, 'YScale', 'log')

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
boxplot2(tmp, 'Width', .6, 'Colors', cmp, 'Labels', ageGroup)
title('Distance')
ylabel('um')
set(findobj(gca,'type','line'),'linew',3)
set(gca,'linew',3, 'FontSize', 22)
% ylim([50000, 1000000])
saveas(h, 'groupedCenterDist_test.png')
saveas(h, 'groupedCenterDist_test.eps')






h = figure;
set(h, 'Position', [0 0 600, 500])
for i = 1:size(ageGroup, 2)
    tmp{i} = squeeze(median_center_stepSpeed(i, :, :));
    tmp{i} = tmp{i}(tmp{i} > 0);
    scatter(i*ones(size(tmp{i})),tmp{i}, [], cmp(i, :), 'LineWidth', 2); hold on
end
set(gca, 'YScale', 'log')

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),'whisker', 2, varargin{:});
boxplot2(tmp, 'Width', .6, 'Colors', cmp, 'Labels', ageGroup)
box off
title('step speed')
ylabel('um')
set(findobj(gca,'type','line'),'linew',3)
set(gca,'linew',3, 'FontSize', 22)
% ylim([50000, 1000000])
saveas(h, 'groupedCenterStepSpeed_test.png')
saveas(h, 'groupedCenterStepSpeed_test.eps')




h = figure;
set(h, 'Position', [0 0 600, 500])
anova_tmp = []; group_tmp = [];
for i = 1:size(ageGroup, 2)
    tmp{i} = squeeze(median_center_stepDist_sum(i, :, :));
    tmp{i} = tmp{i}(tmp{i} > 0);
    anova_tmp = [anova_tmp; tmp{i}];
    group_tmp = [group_tmp; i * ones(length(tmp{i}), 1)];
    scatter(i*ones(size(tmp{i})),tmp{i}, [], cmp(i, :), 'LineWidth', 2); hold on
end
% set(gca, 'YScale', 'log')

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),'whisker', 2, varargin{:});
boxplot2(tmp, 'Width', .6, 'Colors', cmp, 'Labels', ageGroup)
box off
title('step distance sum')
ylabel('um')
set(findobj(gca,'type','line'),'linew',3)
set(gca,'linew',3, 'FontSize', 22)
ylim([0, 2500])
saveas(h, 'groupedCenterStepDistSum_test.png')
saveas(h, 'groupedCenterStepDistSum_test.eps')

% star_distance = zeros(1, size(comb, 1));
% for c = 1 : size(comb, 1)
%     [h_distance(c), p_distance(c)] = ttest2(tmp{comb(c, 1)}, tmp{comb(c, 2)});
% end
% star_distance(p_distance < 0.05) = 1;
% star_distance(p_distance < 0.01) = 2;
% star_distance(p_distance < 0.001) = 3;
% 
[dist_anova_p,dist_t,dist_stats] = anova1(anova_tmp, group_tmp, 'off');
[dist_c, m] = multcompare(dist_stats);
dist_t
for i = 1 : length(tmp)
    aa(1, i) = mean(tmp{i});
    aa(2, i) = std(tmp{i})/sqrt(length(tmp{i}));
end



h = figure;
set(h, 'Position', [0 0 600, 500])
anova_tmp = []; group_tmp = [];
for i = 1:size(ageGroup, 2)
    tmp{i} = squeeze(medianDuration(i, :, :));
    tmp{i} = tmp{i}(tmp{i} > 0);
    anova_tmp = [anova_tmp; tmp{i}];
    group_tmp = [group_tmp; i * ones(length(tmp{i}), 1)];
    scatter(i*ones(size(tmp{i})),tmp{i}, [], cmp(i, :), 'LineWidth', 2); hold on
end
% set(gca, 'YScale', 'log')

col=@(x)reshape(x,numel(x),1);
% boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)), 'whisker', 2, varargin{:});
boxplot2(tmp, 'Width', .6, 'Colors', cmp, 'Labels', ageGroup)
box off
title('Duration')
ylabel('s')
set(findobj(gca,'type','line'),'linew',3)
set(gca,'linew',3, 'FontSize', 22)
ylim([0, 15])
saveas(h, 'groupedDuration_test.png')
saveas(h, 'groupedDuration_test.eps')

% star_duration = zeros(1, size(comb, 1));
% for c = 1 : size(comb, 1)
%     [h_duration(c), p_duration(c)] = ttest2(tmp{comb(c, 1)}, tmp{comb(c, 2)});
% end
% star_duration(p_duration < 0.05) = 1;
% star_duration(p_duration < 0.01) = 2;
% star_duration(p_duration < 0.001) = 3;
[duration_anova_p,duration_t,duration_stats] = anova1(anova_tmp, group_tmp, 'off');
[duration_c, m] = multcompare(duration_stats);
duration_t
for i = 1 : length(tmp)
    aa(1, i) = mean(tmp{i});
    aa(2, i) = std(tmp{i})/sqrt(length(tmp{i}));
end



% plot group interval
h = figure;
set(h, 'Position', [0 0 600, 500])
anova_tmp = []; group_tmp = [];
for i = 1:size(ageGroup, 2)
    tmp{i} = squeeze(medianInterval(i, :, :));
    tmp{i} = tmp{i}(tmp{i} > 0);
    anova_tmp = [anova_tmp; tmp{i}];
    group_tmp = [group_tmp; i * ones(length(tmp{i}), 1)];
    scatter(i*ones(size(tmp{i})),tmp{i}, [], cmp(i, :), 'LineWidth', 2); hold on
end
% set(gca, 'YScale', 'log')

col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
boxplot2(tmp, 'Width', .6, 'Colors', cmp, 'Labels', ageGroup)
box off
title('Interval')
ylabel('s')
set(findobj(gca,'type','line'),'linew',3)
set(gca,'linew',3, 'FontSize', 22)
ylim([0, 160])
saveas(h, 'groupedInterval_test.png')
saveas(h, 'groupedInterval_test.eps')

% star_interval = zeros(1, size(comb, 1));
% for c = 1 : size(comb, 1)
%     [h_interval(c), p_interval(c)] = ttest2(tmp{comb(c, 1)}, tmp{comb(c, 2)});
% end
% star_interval(p_interval < 0.05) = 1;
% star_interval(p_interval < 0.01) = 2;
% star_interval(p_interval < 0.001) = 3;

[interval_anova_p,interval_t,interval_stats] = anova1(anova_tmp, group_tmp, 'off');
[interval_c, m] = multcompare(interval_stats);
interval_t
for i = 1 : length(tmp)
    aa(1, i) = mean(tmp{i});
    aa(2, i) = std(tmp{i})/sqrt(length(tmp{i}));
end




save('data4jointPlot_tra2b.mat', 'medianInterval', 'medianDuration', 'median_center_stepDist_sum', 'medianArea_max', 'individual_var')

save('AnovaPvalue_totalCircularVar.mat', 'area_c', 'dist_c', 'duration_c', 'interval_c', 'total_var', 'total_circ_mean', ...
    'area_t', 'dist_t', 'duration_t', 'interval_t', 'area_anova_p', 'dist_anova_p', 'duration_anova_p', 'interval_anova_p')


save('stats_results.mat', 'star_interval', 'h_interval', 'p_interval', 'star_duration', 'h_duration', 'p_duration', ...
    'star_distance', 'star_distance', 'h_distance', 'star_area', 'star_area', 'h_area', 'area_c', 'dist_c', 'duration_c', 'interval_c')