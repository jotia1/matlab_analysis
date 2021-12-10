function [Suite2p_traces, ROI_centroids] = load_single_fish(pipeline_output_path, fish_number, load_s2p, load_rois)
%% LOAD_SINGLE_FISH - Load s2p and/or ANTs roi results for a fish into matlab
%   Open folders in pipeline_output_path corresponding to fish_number with
%   'suite2p_'/'ants_' and loading relevant suite2p/ants data. Assumes
%   pipeline_output_path contains individual fish folders which in turn
%   contain multiple planeX folders with corresponding suite2p output or
%   or ANTs warped ROI files (ROIs_zbrainspace_XX.csv) 
%
%   Optionally, can specify to load just s2p or just ROI data by setting
%   load_s2p or load_rois to false. Default to true for both. 
%   
%   Args:
%       pipeline_output_path - A full path to a folder containing s2p
%           processed fish.
%       fish_number - zero-padded fish number of the fish to load (e.g. 05)
%       load_s2p - true (default) or false if suite2p data should be loaded
%       load_rois - true (default) or false if ROIs data should be loaded
%   
%
%   Example usage:
%       [Suite2p_traces, ROI_centroids] = load_single_fish('I:\MECP2GEN-Q4070\SPIM\PipelineOutputs', '05');
%   For just suite2p traces:
%       [Suite2p_traces, ROI_centroids] = load_single_fish('I:\MECP2GEN-Q4070\SPIM\PipelineOutputs', '05', true, false);

if ~exist('load_s2p', 'var')
    load_s2p = true;
end
if ~exist('load_rois', 'var')
    load_rois = true;
end

%% Load Suite2p data
fish_folder = dir([pipeline_output_path, '\suite2p_*fish', fish_number, '*']);
Suite2p_traces = []; 
ROI_centroids = [];

% Find number of planes
nplanes = numel(dir([pipeline_output_path, '\suite2p_*fish', fish_number, '*\plane*']));

if load_s2p
    for plane_idx = 0:nplanes-1 % planes index from 0 so subtract 1 from nplanes
        filename = strcat(pipeline_output_path, '\', fish_folder.name, '\plane', num2str(plane_idx), '\Fall.mat'); %access the matfile of data
        load(filename, 'F', 'iscell'); % just load fluorescene traces and centroid data from the matfile

        Suite2p_traces = vertcat(Suite2p_traces, F(iscell(:,1) == 1,:)); % Add the traces to the list
    end
end


if load_rois
    %% Load ants ROIs for this fish
    ants_folder = dir([pipeline_output_path, '\ants_*fish', fish_number, '*']);
    ants_filename = strcat(pipeline_output_path, '\', ants_folder.name, '\ROIs_zbrainspace_', fish_number, '.csv');
    zbrain_rois = readmatrix(ants_filename);
    zbrain_rois(isnan(zbrain_rois)) = 0;
    ROI_centroids = zbrain_rois(:, 1:3);
end

end