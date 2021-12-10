% if sample data is already loaded, dont bother again
if ~exist('spont_df', 'var')
    load('sample_data2.mat');
end

% reorganise ANTs ROIs
ROIs_rounded = zeros(size(ROI_centroids));
ROIs_rounded(:,1) = round(ROI_centroids(:,2));
ROIs_rounded(:,2) = round(ROI_centroids(:,1));
ROIs_rounded(:,3) = round(ROI_centroids(:,3)/2);

% Create plot of rois that are in the fish's brain coloured by region
% You will have to load these masks from your local computer
load('I:\PIPEDATA-Q4414\Zbrain_Masks.mat', 'Zbrain_Masks');
PerBrainRegions = getPerBrainRegions(Zbrain_Masks, ROIs_rounded);
RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Stratum'};
figure
hold on
for i = 1 : length(RegionList)
    region_name = RegionList{i};
    in_region = PerBrainRegions.(region_name).idx;
    plot3(ROI_centroids(in_region, 1), ROI_centroids(in_region, 3), ROI_centroids(in_region, 2), '.')
end





function [ PerBrainRegions ] = getPerBrainRegions(Zbrain_Masks, ROI_correct)

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

end