%% Making a brain Mask 
%load 'Zbrain_Masks'
RegionList_border={'Mesencephalon','Telencephalon', 'Hindbrain','Semicircularis','MON', 'Tectum','Habenula'};
All_Mask=[]
Border_masks=struct;
for i=1:length(RegionList_border);  
    regionName=RegionList_border{i}%;
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
    elseif strcmp(regionName,'Msystem')
        Mask=[];
        Mask_temp=[];
        Msystem_masks=[184 186 187];
        for j=1:3
           %Mask_temp=Zbrain_Masks{Msystem_masks(j),3};
           Mask=vertcat(Mask,Zbrain_Masks{Msystem_masks(j),3});
        end
    
        clear Mask_temp
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
    All_Mask=vertcat(All_Mask,Mask);
    Border_masks.(RegionList_border{i})=Mask;
end

All_Mask_all=unique(All_Mask,'rows');
Brain_Mask= boundary(All_Mask_all(:,1),All_Mask_all(:,2),1); %the 0.8 at the end is the kinda tightness (i guess you could call it) to the plots, a lower value is looser while 1 is tighter
figure;plot(All_Mask_all(Brain_Mask,1),All_Mask_all(Brain_Mask,2))

clear All_Mask Border_masks Hindbrain_Mask IndexC IsInEyes_temp j i Mask RegionList_border regionName Zbrain_Masks 




