function [] = saveas_Raw_data(input_data)
%% Convert pipeline output to variables for Pipeline_analysis_master.m
%   Adjust the variable names to make them as close as possible to
%   variables compatible with the Pipeline_analysis_master.m script
%
%   Inputs:
%       input_data - Path to a file containing raw data loaded from the new
%           pipeline with variables Suite2p_traces, ROI_centroids, 
%           fish_ncells, fish_numbers
%       
%
%
%   Variables in Raw_data.mat:
%       - *File_list 
%       - *File_names 
%       - FishList
%       - ROI_centroids
%       - *ROI_fishnums
%       - ROIs_perfish
%       - *Suite2p_cellNs
%       - Suite2p_traces
%       - ^ZS
%       - *fish_nums
%       - *i
%       - idx_fish
%       - ^traces_detrended
%
%   * indicates the variable is specific to the old method (e.g. filenames
%   of the Slice folders) and cannot be reproduced or does not make sense in
%   the new format. 
%   ^ indicates Maya thought this variable was better replaced by DeltaF
%

load(input_data);

FishList = cellfun(@(x) str2num(x), fish_numbers);
ROIs_perfish = fish_ncells';

idx_fish = zeros(size(Suite2p_traces, 1), 1);
idx = 1;
for fish_idx = 1 : numel(FishList)
    idx_fish(idx: idx + ROIs_perfish(fish_idx) - 1) = zeros([ROIs_perfish(fish_idx), 1]) + FishList(fish_idx);
    idx = idx + ROIs_perfish(fish_idx);
end

save('Raw_data.mat', 'FishList', 'ROI_centroids', 'ROIs_perfish', 'Suite2p_traces', 'idx_fish', '-v7.3')

end