function [] = pipeline_core_analysis(pipeline_output_path, fish_number, sep_idxs)
%% PIPELINE_CORE_ANALYSIS - default analysis and visualisation of pipelin
%   The default per fish analysis and visualisation to perform at the end
%   of the pipeline. 
%
%   Args:
%       pipeline_output_path - Path to output folder of the pipeline
%       fish_number - (string) zero-padded fish number of the fish to load (e.g. 05)
%       sep_idxs - should be the end frame of each stimuli train
%
%   Config structure must contains the following attributes:
%       sep_idxs - list of indices at which stimuli are separated (e.g.
%           spontaneous, auditory, visual trains)
%       regressor - the regressor to use on the auditory stimuli train
%
%   Example usage:
%       pipeline_core_analysis('I:\SCN1LABSYN-Q3714\SPIM\pipeline', '04', [1200]);
%

analysis_dir = fullfile(pipeline_output_path, sprintf('analysis_%s', fish_number));
if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

% if exists, load mat file directly, else do fresh load
save_path = fullfile(analysis_dir, sprintf('raw_fish_%s.mat', fish_number));
if exist(save_path, 'file')
    fprintf('Found existing matlab file, loading that.\n');
    load(save_path);
else
    [Suite2p_traces, stim_trains, ROI_centroids, fish_number] = load_fish_standard_format(pipeline_output_path, fish_number, sep_idxs);
    % Save result so we don't have to recompute
    save(save_path, 'Suite2p_traces', 'stim_trains', 'ROI_centroids', 'fish_number', '-v7.3');
end

%% 
% plot average raw trace and average delta trace
figure('Position', [1, 1, 1920, 1080]);
last_idx = 1;
sep_idxs = [sep_idxs, size(Suite2p_traces, 2)]; % Last idx is seperation
for st = 1 : numel(stim_trains)
    sep_idx = sep_idxs(st);
    subplot(numel(stim_trains), 2, st * 2 - 1);
    plot(mean(Suite2p_traces(:, last_idx:sep_idx),1));
    title(sprintf('Mean raw trace stim train %d (fish%s)', st, fish_number)); ylabel('f');xlabel('frames');
    subplot(numel(stim_trains), 2, st * 2);
    plot(mean(stim_trains{st},1));
    title(sprintf('Mean Df stim train %d (fish%s)', st, fish_number));ylabel('df/f');xlabel('frames');
    last_idx = last_idx + sep_idx;
end
saveas(gcf, fullfile(analysis_dir, sprintf('raw_vs_DF_traces_fish%s.png', fish_number)));

%% Assign each ROI to one of the 11 brain regions
% Create full brain mask
if ismac
    load('/Volumes/UQ-Inst-Gateway1/PIPEDATA-Q4414/Zbrain_Masks.mat');
elseif isunix  % if unix but not mac (Jess) -> should be linux (hpc) or Sarah
    load('/QRISdata/Q4414/Zbrain_Masks.mat')
else  % everyone else
    load('I:\PIPEDATA-Q4414\Zbrain_Masks.mat');
end
PerBrainRegions = getPerBrainRegions(Zbrain_Masks, ROI_centroids);
RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain'};
h = figure;
h.Position = [100, 100, 1180, 1700]./2;
hold on
for i = 1 : length(RegionList)
    region_name = RegionList{i};
    in_region = PerBrainRegions.(region_name).idx;
    plot3(ROI_centroids(in_region, 2), ROI_centroids(in_region, 1), ROI_centroids(in_region, 3), '.')
end
title(sprintf('All ROIs in 11 main regions (fish%s)', fish_number));
saveas(gcf, fullfile(analysis_dir, sprintf('zbrain_warped_ROIs_regional_matlabfig_fish%s.fig', fish_number)));
saveas(gcf, fullfile(analysis_dir, sprintf('zbrain_warped_ROIs_regional_matlabfig_fish%s.png', fish_number)));

%% Draw 2D boundry and take still shots
[tx, ty, sx, sy] = getZbrainOutline(RegionList);
% Add outlines to previous regionwise plot
%plot3(tx, ty, mean(ROI_centroids(:, 3)) * ones(size(ty, 1), 1), 'linewidth', 10); % horizontal

figure; % Front on of fish
plot(ROI_centroids(:, 2), ROI_centroids(:, 3), '.');
hold on
plot(sx, sy)
title(sprintf('All ROIs from fish%s', fish_number));
saveas(gcf, fullfile(analysis_dir, sprintf('zbrain_warped_ROIs_front_fish%s.png', fish_number)));

h = figure; % Top view
h.Position = [100, 100, 1080, 1600]./2;
plot(ROI_centroids(:, 2), ROI_centroids(:, 1), '.');
hold on
plot(tx, ty);
title(sprintf('All ROIs from fish%s', fish_number));
saveas(gcf, fullfile(analysis_dir, sprintf('zbrain_warped_ROIs_top_fish%s.png', fish_number)));

%% for each region, plot the mean
handles = plot_region_wise_means(Zbrain_Masks, ROI_centroids, RegionList, stim_trains);
for h_idx = 1 : numel(handles)
    saveas(handles(h_idx), fullfile(analysis_dir, sprintf('regional_mean_df_stim%d_fish%s.png', h_idx, fish_number)));
end

%% Number of ROIs per region
figure;
n_rois = zeros(numel(RegionList), 1);
total_regional_rois = 0;
for i = 1 : numel(RegionList)
    region_name = RegionList{i};
    n_rois(i) = numel(PerBrainRegions.(region_name).idx);
    total_regional_rois = total_regional_rois + n_rois(i);
end
bar(n_rois);
title(sprintf('ROIs per brain region fish%s (%d total regional ROIs)', fish_number, total_regional_rois));
ax=gca; ax.YAxis.Exponent=0;
xtickangle(45)
xticklabels(RegionList)
xlabel('Region');ylabel('Number of ROIs');
saveas(gca, fullfile(analysis_dir, sprintf('ROIs_per_region_fish%s.png', fish_number)));


%% Regression
regressor = ASD_standard_regressor();

% Plot regressor so we can visually check
figure;
imagesc(regressor);
title(sprintf('regressor used fish%s', fish_number));
saveas(gcf, fullfile(analysis_dir, sprintf('regressor_used_fish%s.png', fish_number)));

% Hardcode regression for stim_train 2 (audio), maybe make this
% configurable later with the config file.
reg_train = stim_trains{2};

% Do Regression 
model_basic=struct();
parfor i=1:size(reg_train,1)
    mdl=fitlm(regressor',reg_train(i,:));  %change here to use 2 regressor Hab(:,220 bla bla bla) 
    model_basic(i).coef=mdl.Coefficients;        
    model_basic(i).rsquared=mdl.Rsquared.Adjusted;
end

%See what the r2 looks like
rsq=[model_basic.rsquared];
figure;
histogram(rsq);
title(sprintf('All r2 values for linear regression fish%s', fish_number));
xlabel('r^2'); ylabel('count'); 
saveas(gcf, fullfile(analysis_dir, sprintf('r2_values_fish%s.png', fish_number)));


%% Now some regression analysis

idx_rsq=find(rsq>0.05);

% Look at where the neurons are that pass the rsq filter
figure('Position', [100, 100, 1080, 1600]./2);
hold on
for region_idx = 1 : numel(RegionList)
    region_name = RegionList{region_idx};
    region_rois = PerBrainRegions.(region_name).idx;
    shared = intersect(idx_rsq, region_rois);
    plot3(ROI_centroids(shared,2),ROI_centroids(shared,1),ROI_centroids(shared,3),'.');
end
% Draw outline
plot(tx, ty);
saveas(gcf, fullfile(analysis_dir, sprintf('aud_responsive_roi_locs_fish%s.png', fish_number)));

% Mean of all auditory responsive neurons
figure;
plot(mean(reg_train(idx_rsq, :)));
title(sprintf('Mean of auditory responsive neurons, fish%s', fish_number));
xlabel('Frame');
ylabel('mean df/f');
saveas(gcf, fullfile(analysis_dir, sprintf('aud_responsive_mean_df_fish%s.png', fish_number)));

% Mean of auditory responsive neurons by region
figure('Position', [1, 1, 1920, 1080]);
for region_idx = 1 : numel(RegionList)
    subplot(4, 3, region_idx);
    region_name = RegionList{region_idx};
    region_rois = PerBrainRegions.(region_name).idx;
    shared = intersect(idx_rsq, region_rois);
    region_resp_roi_dfs = reg_train(shared, :);
    plot(mean(region_resp_roi_dfs));
    title(sprintf('%s mean df/f', region_name));
    xlabel('Frame'); ylabel('df/f');
end
sgtitle(sprintf('Region wise mean of audio responsive rois fish%s', fish_number));
saveas(gcf, fullfile(analysis_dir, sprintf('aud_responsive_regional_mean_df_fish%s.png', fish_number)));

