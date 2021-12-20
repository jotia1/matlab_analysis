%% Not meant as a feature of this code, individual fish should already exist 
% when pipeline does per fish analysis, early versions of the code were
% missing this functionality so I added this hack in here to come back and
% add the individual fish matlab files posthoc. 211220 - Josh


pipeline_output_path = 'I:\SCN1LABSYN-Q3714\SPIM\pipeline';
sep_idxs = [1200];

fish_folders = dir([pipeline_output_path, '\suite2p_*']);

num_fish = numel(fish_folders);

%% Get all fish numbers, padded with leading zeros (e.g. 05 rather than 5)
fish_folder_names = {fish_folders.name};
fin = cellfun(@(x)regexp(x,'fish(\d+)','tokens'), fish_folder_names, 'UniformOutput', false);
fish_numbers = cell(numel(fin), 1);
for i = 1:numel(fin)
    fish_numbers{i} = fin{i}{1}{1};
end


for fish_idx = 1:num_fish
    fish_number = fish_numbers{fish_idx};
    matfile_name = fullfile(pipeline_output_path, sprintf('analysis_%s', fish_number), sprintf('raw_fish_%s.mat', fish_number));
    if exist(matfile_name, 'file')
        fprintf('Found existing matlab file, skipping (fish%s)\n', fish_number)
        continue
    end
    
    pipeline_core_analysis(pipeline_output_path, fish_number, sep_idxs);
    
end

