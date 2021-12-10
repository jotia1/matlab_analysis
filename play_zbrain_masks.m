%% 

if ~exist('Zbrain_Masks', 'var')
    load('I:\PIPEDATA-Q4414\Zbrain_Masks.mat');
end


for i = 1:size(Zbrain_Masks, 1)
    region = Zbrain_Masks{i, 3};
    
end

