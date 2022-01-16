function [ROI_centroids, fish_ncells, fish_numbers, stim_trains] = load_all_fish_standard_format(pipeline_output_path, sep_idxs)
%% LOAD_ALL_FISH_STANDARD_FORMAT - Load all fish s2p and/or ANTs roi results into matlab
%
%
%
%
%   Example usage:
%       


% set up default values (true) for load_s2p/rois, raise error if both false
if ~exist('load_s2p', 'var')
    load_s2p = true;
end
if ~exist('load_ants', 'var')
    load_rois = true;
end
if ~load_s2p && ~load_ants
    throw(MException('LOAD_ALL_FISH:NothingToLoad', 'load_s2p and load_rois cannot both be false.'))
end

fish_folders = dir([pipeline_output_path, '\suite2p_*']);

% TODO : testing hack
%fish_folders = fish_folders(1:2);
%fish_folders(22) = [];  % Remove fish 41 which is missing data


num_fish = numel(fish_folders);

%% Get all fish numbers, padded with leading zeros (e.g. 05 rather than 5)
fish_folder_names = {fish_folders.name};
fin = cellfun(@(x)regexp(x,'fish(\d+)','tokens'), fish_folder_names, 'UniformOutput', false);
fish_numbers = cell(numel(fin), 1);
for i = 1:numel(fin)
    fish_numbers{i} = fin{i}{1}{1};
end


%% Loop through folders to get traces and xy locations of all ROIs Suite2p defines as cells
%Suite2p_traces = []; 
ROI_centroids = [];
stim_trains = cell(numel(sep_idxs) + 1, 1);
fish_ncells = zeros(num_fish, 1); % number of cells per fish

for fish_idx = 1:num_fish
    folder = fish_folders(fish_idx).name;
    
    % fish41 is missing data from mecp2
    if strcmp(folder(19:end), 'fish41') == 1
        continue  % Skip fish
    end
    
    fish_number = fish_numbers{fish_idx};
    matfile_name = fullfile(pipeline_output_path, sprintf('analysis_%s', fish_number), sprintf('raw_fish_%s.mat', fish_number));
    if exist(matfile_name, 'file')
        fprintf('Found existing matlab file, loading that (fish%s)\n', fish_number)
        data = load(matfile_name,'stim_trains', 'ROI_centroids');
        fish_stim_trains = data.stim_trains;
        fish_ROI_centroids = data.ROI_centroids;
    else
        [~, fish_stim_trains, fish_ROI_centroids, ~] = load_fish_standard_format(pipeline_output_path, fish_number, sep_idxs);
    end
    
    % Join individual fish with the collective fish
    %Suite2p_traces = vertcat(Suite2p_traces, fish_Suite2p_traces);
    ROI_centroids = vertcat(ROI_centroids, fish_ROI_centroids);
    ROI_centroids(isnan(ROI_centroids)) = 0; % Avoid that nan ROI at end of file from ANTs
    for st = 1 : numel(stim_trains)
        stim_trains{st} = vertcat(stim_trains{st}, fish_stim_trains{st});
    end
    
    % Count cells as #traces if loading traces, else use #ROIs
    ncells = size(fish_stim_trains{1}, 1);
    if ~load_s2p
        ncells = size(fish_ROI_centroids, 1);
    end
    fish_ncells(fish_idx) = ncells;
    
end



end