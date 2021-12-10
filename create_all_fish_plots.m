function [] = create_all_fish_plots(pipeline_output_path)

[Suite2p_traces, ROI_centroids, fish_ncells, fish_numbers] = load_all_fish(pipeline_output_path);

%% Assume all data will be stitched, seperate into spont and aud
% Remove linear trend from data to account for decrease in baseline fluorescence
spont_traces_detrended = detrend(Suite2p_traces(1:600, :), 'linear')';
aud_traces_detrended = detrend(Suite2p_traces(601:end, :), 'linear')';

% Z-score the detrended data
spont_ZS = zscore(spont_traces_detrended, 1, 2);
aud_ZS = zscore(aud_traces_detrended, 1, 2);


%% Create plots of raw, detrended and zscored traces
figure; 
subplot(2,2,1); plot(mean(Suite2p_traces(1:600),1)); title('Mean of raw traces (spont)'); ylabel('df/f');xlabel('frames')
subplot(2,2,2); plot(mean(spont_traces_detrended,1)); title('Mean of detrended traces (spont)');ylabel('df/f (after detrending)');xlabel('frames')
subplot(2,2,3); plot(mean(spont_ZS,1)); title('Mean of all Z-scored traces (spont)'); ylabel('z-scored df/f');xlabel('frames')
%saveas(gcf,'spont_raw_data_all_ROIs.png')

figure; 
subplot(2,2,1); plot(mean(Suite2p_traces(601:end),1)); title('Mean of raw traces (aud)'); ylabel('df/f');xlabel('frames')
subplot(2,2,2); plot(mean(aud_traces_detrended,1)); title('Mean of detrended traces (aud)');ylabel('df/f (after detrending)');xlabel('frames')
subplot(2,2,3); plot(mean(aud_ZS,1)); title('Mean of all Z-scored traces (aud)'); ylabel('z-scored df/f');xlabel('frames')
%saveas(gcf,'aud_raw_data_all_ROIs.png')

figure; 
bar(fish_ncells); xticklabels(fish_numbers); xlabel('Fish number'); ylabel('Number of ROIs'); 
title('Total number of ROIs per fish'); set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')));
%saveas(gcf,'Total_ROIs_per_fish.png');


%% Do zbrains stuff
load('I:\PIPEDATA-Q4414\MaskDatabase.mat');

% reorganise ANTs ROIs
ROIs = zeros(size(ROI_centroids));
ROIs(:,1) = height - round(ROI_centroids(:,2));
ROIs(:,2) = round(ROI_centroids(:,1));
ROIs(:,3) = round(ROI_centroids(:,3)/2);

% Create full brain mask
full_brain_mask = zeros(size(MaskDatabase, 1), 1);
region_wise_masks = cell(294, 1);
for region_idx = [1:77 79:294]  % Exlcude eyes
    full_brain_mask = full_brain_mask | MaskDatabase(:, region_idx);
    linear_locations = find(MaskDatabase(:, region_idx));
    [rx, ry, rz] = ind2sub([height, width, Zs], linear_locations);
    region_wise_masks{region_idx} = [rx, ry, rz];
end
% TODO : could replace the below with a vertcat of region_wise_masks (some
% unnecessary computation here.)
f = find(full_brain_mask);
[mx, my, mz] = ind2sub([height, width, Zs], f);
mask_ROIs = [mx, my, mz]; % Now full brain mask is an [X, 3] matrix of coords

% select only ROIs that are in the brain
in_brain = ismember(ROIs, mask_ROIs, 'rows');
rois_inbrain = ROIs(in_brain, :);

% Create plot of ALL in brain ROIs from ALL fish
plot3(rois_inbrain(:, 1), rois_inbrain(:, 2), rois_inbrain(:, 3), '.')
title('All ROIs from all fish');
% Save a matlab figure so we can interact with plot later
saveas(gcf, 'all_ROI_locations_all_fish.fig');

% Now save some png side views
title('All ROIs from all fish top');
view(0, 90);
%saveas(gcf, 'all_ROI_locations_all_fish_top.png');
title('All ROIs from all fish front');
view(270, 0);
%saveas(gcf, 'all_ROI_locations_all_fish_front.png');
title('All ROIs from all fish side');
view(0, 0);
%saveas(gcf, 'all_ROI_locations_all_fish_side.png');


% Sort ROIs per region
PerBrainRegions=struct();
RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','Hindbrain','Stratum'};
for i=1:length(RegionList)
    regionName=RegionList{i};
    if strcmp(regionName,'Telencephalon')
        Mask = region_wise_masks{294};
    elseif strcmp(regionName,'Hindbrain')
        Hindbrain_Mask=region_wise_masks{259};
        Mask=region_wise_masks{131};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        %Mask=region_wise_masks{295};
        %IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
        %Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Hindbrain_Mask;
    else
        Mask=[];
        IndexC=strfind({MaskDatabaseNames{1, :}}, regionName);
        IndexC=find(not(cellfun('isempty', IndexC)));
        for j=IndexC
            if isempty(Mask)
                Mask=region_wise_masks{j};
            else
                Mask=vertcat(Mask,region_wise_masks{j});
            end
        end
    end
    Mask=unique(Mask,'rows');
    IsInBrainRegion=ismember(rois_inbrain,Mask,'rows');
    PerBrainRegions.(regionName).idx=find(IsInBrainRegion==1);       
end

ROIs_perfish = fish_ncells';
idx_Fish = zeros(size(Suite2p_traces, 1), 1);
idx = 1;
for fish_idx = 1 : numel(FishList)
    idx_Fish(idx: idx + ROIs_perfish(fish_idx) - 1) = zeros([ROIs_perfish(fish_idx), 1]) + FishList(fish_idx);
    idx = idx + ROIs_perfish(fish_idx);
end

% Now plot histograms of regions
progressbar;
FishList = cellfun(@(x) str2num(x), fish_numbers);
h=zeros(length(RegionList),length(FishList));
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
    regionName=RegionList{i};
    idx_temp=PerBrainRegions.(regionName).idx;    
    for fish=1:length(FishList)
        temp_WT_fish=length(find(idx_Fish(idx_temp)==FishList(fish)));
        h(i,fish)=temp_WT_fish;
    end    
end
figure;bar(sum(h,2));




end