function [] = create_single_fish_plots(pipeline_output_path, fish_number)
%% 
%
%
%       fish_number should be zero-padded (e.g. 05 rather than 5)

[Suite2p_traces, ROI_centroids] = load_single_fish(pipeline_output_path, fish_number);


%% Assume all data will be stitched, seperate into spont and aud
sep_idx = 1200; % 600 sec * 2 fps

% Delta traces
spont_df = DeltaF2(Suite2p_traces(:, 1:sep_idx),101,7);
spont_df(isnan(spont_df))=0;
aud_df = DeltaF2(Suite2p_traces(:, sep_idx+1:end),101,7);
aud_df(isnan(aud_df))=0;

figure; % Average raw trace and average delta trace
subplot(2,2,1); plot(mean(Suite2p_traces(:, 1:sep_idx),1)); title(sprintf('Mean raw trace spont (fish%s)', fish_number)); ylabel('f');xlabel('frames')
subplot(2,2,2); plot(mean(spont_df,1)); title(sprintf('Mean Df spont (fish%s)', fish_number));ylabel('df/f');xlabel('frames')
subplot(2,2,3); plot(mean(Suite2p_traces(:, sep_idx+1:end),1)); title(sprintf('Mean raw trace aud (fish%s)', fish_number)); ylabel('f');xlabel('frames')
subplot(2,2,4); plot(mean(aud_df,1)); title(sprintf('Mean Df aud (fish%s)', fish_number));ylabel('df/f');xlabel('frames')
saveas(gcf,sprintf('raw_vs_DF_traces_fish%s.png', fish_number))

%% Do zbrains stuff
load('I:\PIPEDATA-Q4414\MaskDatabase.mat', 'height', 'MaskDatabase', 'MaskDatabaseNames', 'width', 'Zs');

% reorganise ANTs ROIs
ROIs_rounded = zeros(size(ROI_centroids));
ROIs_rounded(:,1) = round(ROI_centroids(:,2));
ROIs_rounded(:,2) = round(ROI_centroids(:,1));
ROIs_rounded(:,3) = round(ROI_centroids(:,3)/2);


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


% Create plot of rois that are in the fish's brain coloured by region
h = figure;
h.Position = [100, 100, 1920, 1080]./2;
hold on
ROI_regions = zeros(size(ROI_centroids, 1), 1);
for region_idx = [1:77 79:294]  % Exlcude eyes
    in_region = ismember(ROIs_rounded, region_wise_masks{region_idx}, 'rows');
    ROI_regions(in_region) = region_idx;
    plot3(ROI_centroids(in_region, 1), ROI_centroids(in_region, 3), ROI_centroids(in_region, 2), '.')
end

title(sprintf('All ROIs from fish%s', fish_number));
% Save a matlab figure so we can interact with plot later
saveas(gcf,sprintf('zbrain_warped_ROIs_matlabfig_fish%s.fig', fish_number));
view(0, 90);
saveas(gcf, sprintf('zbrain_warped_ROIs_front_fish%s.png', fish_number));
view(90, 0);
saveas(gcf, sprintf('zbrain_warped_ROIs_side_fish%s.png', fish_number));
view(180, 0);
h.Position = [100, 100, 1080, 1600]./2;
saveas(gcf, sprintf('zbrain_warped_ROIs_top_fish%s.png', fish_number));


%% Plot region wise mean df traces
load('I:\PIPEDATA-Q4414\Zbrain_Masks.mat');
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3}); %Makes a massive matrix of all the Zbrain regions except the eyes
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

%You can check the alignments with :
idx_rand=randsample(length(Zbrain_AllMask),10000); %need to select a few members only
figure;
subplot(1,2,1);
scatter(Zbrain_AllMask(idx_rand,1),Zbrain_AllMask(idx_rand,2),'.')
hold on;scatter(ROIs_rounded(:,1),ROIs_rounded(:,2),'.');
subplot(1,2,2);
scatter(Zbrain_AllMask(idx_rand,1),Zbrain_AllMask(idx_rand,3),'.')
hold on;scatter(ROIs_rounded(:,1),ROIs_rounded(:,3),'.');
title(sprintf('ANTs alignment check fish%s', fish_number));
saveas(gcf, sprintf('ANTs_alignment_fish%s.png', fish_number));


%% Assign all ROIs to a region
PerBrainRegions=struct();
RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Stratum'};
progressbar;
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
    regionName=RegionList{i};
    if strcmp(regionName,'Telencephalon')
        Mask=Zbrain_Masks{294,3};
        
    elseif strcmp(regionName,'Hindbrain')
        Hindbrain_Mask=Zbrain_Masks{259,3};
        Mask=Zbrain_Masks{131,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Zbrain_Masks{295,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Hindbrain_Mask;
        
    else
        Mask=[];
        IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
        IndexC=find(not(cellfun('isempty', IndexC)));
        for j=IndexC
            if isempty(Mask)
                Mask=Zbrain_Masks{j,3};
            else
                Mask=vertcat(Mask,Zbrain_Masks{j,3});
            end
        end
    end
    Mask=unique(Mask,'rows');
    IsInBrainRegion=ismember(ROI_correct,Mask,'rows'); 
    PerBrainRegions.(regionName).idx=find(IsInBrainRegion==1);    
end


%% for each region, plot the mean
figure;
hold on
for i = 1: numel(RegionList)
    subplot(4, 3, i);
    region_name = RegionList{i};
    region_rois = PerBrainRegions.(region_name).idx;
    region_aud_dfs = aud_df(region_rois, :);
    plot(mean(region_aud_dfs))
    title(sprintf('%s mean df/f', region_name));
    xlabel('Frame'); ylabel('df/f');
end




% Make gif of fish
axis tight manual % this ensures that getframe() returns a consistent size
filename = sprintf('warped_ROIs_fish%s.gif', fish_number);
axis off
title(sprintf('warped ROIs for fish%s', fish_number))

%for n = pi:pi/128:3 *pi
for n = 0:360/256:360
    
    view(n, 0);
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 0
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',1/24); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',1/24); 
    end
end



%% create region-wise mean df/f traces



%% prepare for regression

%Building the regressor by ID the stim times and then applying to data
spike=[0,1.69644104899772,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
spike=spike/max(spike);
Regressor=zeros(1,size(aud_df,2));

%classify timepoints of stimulus - Vol 4 is the loudest Vol 1 is the
%quietist, then there is PPI and Habituation

Vol4=[200, 400, 520, 840, 1080]; %-6dB
for i=Vol4
    Regressor(1,i-1:i-1+length(spike)-1)=spike';
end
Vol3=[120, 440, 680, 920, 1040]; %-12dB
for i=Vol3
    Regressor(1,i-1:i-1+length(spike)-1)=spike';
end
Vol2=[160, 320, 600, 720, 960]; %-18dB
for i=Vol2
    Regressor(1,i-1:i-1+length(spike)-1)=spike';
end
Vol1=[280, 560, 760, 880, 1000]; %-24dB
for i=Vol1
    Regressor(1,i-1:i-1+length(spike)-1)=spike';
end
PPI=[240, 360, 480, 640, 800];
for i=PPI
    Regressor(2,i-1:i-1+length(spike)-1)=spike';
end

AudHab=(10 * 60 * 2 : 6 : 11 * 60 * 2 - 6); % Only includes first habituation block
% second block = (12 * 60 * 2 : 6 : 12.5 * 60 * 2 -6)
for i=1:numel(AudHab)
    Regressor(3,AudHab(i):AudHab(i)+length(spike)-1)=spike';
end

%Regressor=Regressor(:,1:1440);
figure;plot(Regressor')

save('Regressor.mat','Regressor','Vol4','Vol3', 'Vol2', 'Vol1', 'PPI', 'AudHab');

% Do Regression 
model_basic=struct();
parfor i=1:size(aud_df,1)
    mdl=fitlm(Regressor',aud_df(i,:));  %change here to use 2 regressor Hab(:,220 bla bla bla) 
    model_basic(i).coef=mdl.Coefficients;        
    model_basic(i).rsquared=mdl.Rsquared.Adjusted;
end

%See what the r2 looks like
rsq=[model_basic.rsquared];
figure;histogram(rsq);

end