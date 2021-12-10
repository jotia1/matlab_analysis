function [ PerBrainRegions ] = getPerBrainRegions(Zbrain_Masks, ROI_correct)

% Need to make sure ROIs are appropriately rotated to match the zbrain
% masks

% Round just in case the ROIs weren't already rounded
ROI_correct = round(ROI_correct);

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