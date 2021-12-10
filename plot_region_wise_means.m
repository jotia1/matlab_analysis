function [ handles ] = plot_region_wise_means(Zbrain_Masks, ROI_centroids, RegionList, stim_trains)

PerBrainRegions = getPerBrainRegions(Zbrain_Masks, ROI_centroids);
handles = [];
for st = 1 : numel(stim_trains)
    handles(end + 1) = figure('Position', [1, 1, 1920, 1080]);
    hold on
    for i = 1: numel(RegionList)
        subplot(4, 3, i);
        region_name = RegionList{i};
        region_rois = PerBrainRegions.(region_name).idx;
        stim_train = stim_trains{st};
        region_aud_dfs = stim_train(region_rois, :);
        plot(mean(region_aud_dfs))
        title(sprintf('%s mean df/f', region_name));
        xlabel('Frame'); ylabel('df/f');
    end
    sgtitle(sprintf('Stimulus train %d', st));
end


end
