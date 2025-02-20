% generating toy data
rng("default");

% Parameters
n_subjects = 20;
n_rois = 5;
n_timepoints = 50;
n_sessions = 2;

% Generate ROI data (n, r, t, s)
roi_data = randn(n_subjects, n_rois, n_timepoints, n_sessions);  % Normally distributed random data

% Generate behavioral data (1, n) 
behavior = randi([0, 10], 1, n_subjects); % Generate integer scores between 0 and 10
% Specify the data folder path
data_folder = 'data';

% Check if the 'data' folder exists, and create it if it doesn't
if ~isfolder(data_folder)
    mkdir(data_folder);
    disp(['Created folder: ', data_folder]);
end

% Save the data
save(fullfile(data_folder, 'toy_fmri_data.mat'), 'roi_data', 'behavior','n_subjects','n_rois','n_timepoints','n_sessions');

disp('Toy fMRI and behavioral data generated and saved to ./data/toy_fmri_data.mat');

%% %% get_mean_roi_timeseries 
function mean_timeseries = get_mean_roi_timeseries(roi_data)
%GET_MEAN_ROI_TIMESERIES Calculates the mean ROI timeseries for each
%subject and ROI, averaged across sessions.
%
%   mean_timeseries = GET_MEAN_ROI_TIMESERIES(roi_data)
%
%   Input:
%       roi_data: A 4-D (n,r,t,s) matrix where
%           n = number of subjects
%           r = number of ROIs
%           t = number of time points
%           s = number of scan sessions per subject
%
%   Output:
%       mean_timeseries: A 3-D (n,r,t) matrix of mean ROI timeseries

    n_subjects = size(roi_data, 1);
    n_rois = size(roi_data, 2);
    n_timepoints = size(roi_data, 3);

    mean_timeseries = zeros(n_subjects, n_rois, n_timepoints);

    for subj = 1:n_subjects
        for roi = 1:n_rois
            % Average across sessions
            mean_timeseries(subj, roi, :) = mean(roi_data(subj, roi, :, :), 4);
        end
    end

end
%% temporal_isc_pairwise

% Specify the data folder path
data_folder = 'data';

% Load the data
load(fullfile(data_folder, 'toy_fmri_data.mat'), 'roi_data','n_subjects','n_rois','n_timepoints','n_sessions');

% Calculate mean ROI timeseries across sessions
mean_roi_timeseries = get_mean_roi_timeseries(roi_data);

% Preallocate
pairwise_temporal_ISC = zeros(n_subjects, n_subjects, n_rois);  % (n, n, r)


% Temporal ISC calculation
for roi = 1:n_rois
    for subj1 = 1:n_subjects
        for subj2 = 1:n_subjects
            if subj1 ~= subj2
                timeseries1 = squeeze(mean_roi_timeseries(subj1, roi, :));
                timeseries2 = squeeze(mean_roi_timeseries(subj2, roi, :));

                pairwise_temporal_ISC(subj1, subj2, roi) = corr(timeseries1, timeseries2);
            else
                pairwise_temporal_ISC(subj1, subj2, roi) = NaN;
            end
        end
    end
end

% Save the results
save(fullfile(data_folder, 'temporal_isc_pairwise_results.mat'), 'pairwise_temporal_ISC');

disp('Pairwise temporal ISC calculated and saved to ./data/temporal_isc_pairwise_results.mat');
disp(['pairwise_temporal_ISC size: ', num2str(size(pairwise_temporal_ISC))]); %Disp command

%% temporal_isc_leaveoneout

% Specify the data folder path
data_folder = 'data';

% Load the data
load(fullfile(data_folder, 'toy_fmri_data.mat'), 'roi_data','n_subjects','n_rois','n_timepoints','n_sessions');

% Calculate mean ROI timeseries across sessions
mean_roi_timeseries = get_mean_roi_timeseries(roi_data);

% Preallocate
loo_temporal_ISC = zeros(n_subjects, n_rois);  % (n, r)


% Leave-one-out ISC calculation
for roi = 1:n_rois
    for subj_leaveout = 1:n_subjects
        % Create the group timeseries (average of all other subjects)
        group_timeseries = squeeze(mean(mean_roi_timeseries(setdiff(1:n_subjects, subj_leaveout), roi, :), 1));

        % Get the left-out subject's timeseries
        subject_timeseries = squeeze(mean_roi_timeseries(subj_leaveout, roi, :));

        % Calculate the correlation
        loo_temporal_ISC(subj_leaveout, roi) = corr(subject_timeseries, group_timeseries);
    end
end


% Save the results
save(fullfile(data_folder, 'temporal_isc_leaveoneout_results.mat'), 'loo_temporal_ISC');

disp('Leave-one-out temporal ISC calculated and saved to ./data/temporal_isc_leaveoneout_results.mat');
disp(['loo_temporal_ISC size: ', num2str(size(loo_temporal_ISC))]);


%% % dynamic_isc

% Specify the data folder path
data_folder = 'data';

% Load the data
load(fullfile(data_folder,'toy_fmri_data.mat'), 'roi_data','n_subjects','n_rois','n_timepoints','n_sessions');

% Calculate mean ROI timeseries across sessions
mean_roi_timeseries = get_mean_roi_timeseries(roi_data);

window_size = 10;
step_size = 10;
n_windows = floor((n_timepoints - window_size) / step_size) + 1;

% Preallocate
loo_dynamic_ISC = zeros(n_subjects, n_rois, n_windows);  % (n, r, t/w)

% Dynamic ISC calculation
for roi = 1:n_rois
    for subj_leaveout = 1:n_subjects
        for window = 1:n_windows
            start_time = (window - 1) * step_size + 1;
            end_time = start_time + window_size - 1;

            % Extract timeseries for the current window
            subject_timeseries = squeeze(mean_roi_timeseries(subj_leaveout, roi, start_time:end_time));

            % Create group timeseries (average of all other subjects)
            group_timeseries = squeeze(mean(mean_roi_timeseries(setdiff(1:n_subjects, subj_leaveout), roi, start_time:end_time), 1));

            % Calculate correlation
            loo_dynamic_ISC(subj_leaveout, roi, window) = corr(subject_timeseries, group_timeseries);
        end
    end
end


% Save the results
save(fullfile(data_folder,'dynamic_isc_results.mat'), 'loo_dynamic_ISC', 'window_size', 'step_size');

% Plot change in dynamic ISC for the first subject and first ROI
figure;  % Create a new figure
plot(1:n_windows, squeeze(loo_dynamic_ISC(1, 1, :)));
xlabel('Time Window');
ylabel('Dynamic ISC');  % Changed from Temporal ISC to Dynamic ISC
title('Change in Dynamic ISC (Subject 1, ROI 1)');  % Changed from Temporal ISC to Dynamic ISC

% Save the plot to a file
saveas(gcf, fullfile(data_folder,'dynamic_isc_plot.pdf'), 'pdf');

% Show the Mean dynamic ISC value
disp(['Mean Dynamic ISC (across subjects, ROIs, and windows): ' num2str(mean(loo_dynamic_ISC(:)))]);

%% % spatial_isc

% Specify the data folder path
data_folder = 'data';

% Load the data
load(fullfile(data_folder,'toy_fmri_data.mat'), 'roi_data','n_subjects','n_rois','n_timepoints','n_sessions');

% Calculate mean ROI timeseries across sessions
mean_roi_timeseries = get_mean_roi_timeseries(roi_data);

% Preallocate
loo_spatial_ISC = zeros(n_subjects, 1);  % (n, 1)

% Spatial ISC calculation
for subj_leaveout = 1:n_subjects
    % Get mean activity for each ROI, averaged over time, for the left-out subject
    subject_roi_means = squeeze(mean(mean_roi_timeseries(subj_leaveout, :, :), 3));  % (1, r)

    % Get mean activity for each ROI, averaged over time and over all *other* subjects
    group_roi_means = squeeze(mean(mean(mean_roi_timeseries(setdiff(1:n_subjects, subj_leaveout), :, :), 1), 3)); % (1, r)

    % Calculate correlation between the subject's spatial pattern and the group's spatial pattern
    loo_spatial_ISC(subj_leaveout) = corr(subject_roi_means', group_roi_means');
end

% Save the results
save(fullfile(data_folder,'spatial_isc_results.mat'), 'loo_spatial_ISC');

disp('Spatial ISC calculated and saved to ./data');
disp(['loo_spatial_ISC size: ', num2str(size(loo_spatial_ISC))]);

%% intra_subject_correlation

% Specify the data folder path
data_folder = 'data';

% Load the data
load(fullfile(data_folder,'toy_fmri_data.mat'), 'roi_data','n_subjects','n_rois','n_timepoints','n_sessions');  % Make sure the path is correct

n_subjects = size(roi_data, 1);
n_rois = size(roi_data, 2);
n_timepoints = size(roi_data, 3);
n_sessions = size(roi_data, 4);

% Preallocate
intrasubject_temporal_ISC = zeros(n_subjects, n_rois);  % (n, r)

% Intra-subject correlation calculation
for subj = 1:n_subjects
    for roi = 1:n_rois
        % Extract timeseries for session 1 and session 2 for the current subject and ROI
        timeseries_session1 = squeeze(roi_data(subj, roi, :, 1));  % Session 1
        timeseries_session2 = squeeze(roi_data(subj, roi, :, 2));  % Session 2

        % Calculate the correlation between the two sessions
        intrasubject_temporal_ISC(subj, roi) = corr(timeseries_session1, timeseries_session2);
    end
end

% Save the results
save(fullfile(data_folder,'intra_subject_correlation_results.mat'), 'intrasubject_temporal_ISC'); % Make sure the path is correct

disp('Intra-subject correlation calculated and saved to ./data/intra_subject_correlation_results.mat');
disp(['intrasubject_temporal_ISC size: ', num2str(size(intrasubject_temporal_ISC))]);

%% inter_subject_functional_connectivity

% Specify the data folder path
data_folder = 'data';

% Load the data
load(fullfile(data_folder,'toy_fmri_data.mat'), 'roi_data','n_subjects','n_rois','n_timepoints','n_sessions');

% Calculate mean ROI timeseries across sessions
mean_roi_timeseries = get_mean_roi_timeseries(roi_data);

% Preallocate
loo_ISFC = zeros(n_subjects, n_rois, n_rois);  % (n, r, r)

% Inter-subject functional connectivity calculation
for subj_leaveout = 1:n_subjects
    % Calculate mean timeseries for all *other* subjects
    group_timeseries = squeeze(mean(mean_roi_timeseries(setdiff(1:n_subjects, subj_leaveout), :, :), 1));  % (r, t)

    % Calculate ISFC for the left-out subject
    for roi1 = 1:n_rois
        for roi2 = 1:n_rois
            % Extract timeseries for the left-out subject
            subject_timeseries_roi1 = squeeze(mean_roi_timeseries(subj_leaveout, roi1, :));

            % Extract timeseries for the *group*
            group_timeseries_roi2 = squeeze(group_timeseries(roi2, :));  % Keep as row vector

            % Calculate correlation
            loo_ISFC(subj_leaveout, roi1, roi2) = corr(subject_timeseries_roi1, group_timeseries_roi2');
        end
    end
end

% Save the results
save(fullfile(data_folder,'inter_subject_functional_connectivity_results.mat'), 'loo_ISFC'); % Correct path

disp('Inter-subject functional connectivity calculated and saved to ./data/inter_subject_functional_connectivity_results.mat');
disp(['loo_ISFC size: ', num2str(size(loo_ISFC))]);
% Average over ROI pairs and subjects
disp(['Mean ISFC (across subjects and ROI pairs): ' num2str(mean(loo_ISFC(:)))]);

%% representational_similarity_to_behavior 

% In this section, I compute the behavioral similarity matrix and 
% relate it to the temporal inter-subject correlation (ISC) matrix. 
% 
% The behavioral similarity is calculated using the inverse of the 
% absolute difference between each pair of subjects' behavioral scores.
% A higher similarity is assigned to pairs of subjects with closer 
% behavioral scores. This operationalizes the idea that more similar 
% behaviors (i.e., closer behavioral scores) should correspond to 
% more similar neural responses (higher ISC) in relevant brain regions.
%
% Hypothesis:
% I hypothesize that higher behavioral similarity between pairs of 
% subjects (reflected in the behavioral_similarity matrix) will 
% correspond to higher ISC in specific regions of interest (ROIs). 
% This would suggest that behavioral similarity influences the degree 
% of neural synchronization, supporting the idea that ISC may reflect
% shared cognitive or perceptual processes related to behavioral patterns.


% Specify the data folder path
data_folder = 'data';

% Load the data
load(fullfile(data_folder,'toy_fmri_data.mat'), 'roi_data', 'behavior','n_subjects','n_rois','n_timepoints','n_sessions');
load(fullfile(data_folder, 'temporal_isc_pairwise_results.mat'), 'pairwise_temporal_ISC');

% Define shorthand variables
n = n_subjects;
r = n_rois;

% Ensure behavior is (1, n_subjects) or (n_subjects, 1)
if size(behavior, 1) == 1 && size(behavior, 2) == n
    % It's already a row vector, so do nothing
elseif size(behavior, 2) == 1 && size(behavior, 1) == n
    behavior = behavior'; % Transpose to make it a row vector
else
    error('Behavior variable has incorrect dimensions. Expected (1, n_subjects) or (n_subjects, 1).');
end

% 6. Relate inter-subject representational similarity patterns to behavior
% Create a (n,n) matrix for behavioral similarity
behavioral_similarity = zeros(n, n);

% Compute similarity: higher similarity means smaller difference in behavioral scores
for i = 1:n
    for j = 1:n
        behavioral_similarity(i, j) = 1 / (1 + abs(behavior(i) - behavior(j)));
    end
end

% Display result size
disp('Behavioral Similarity Matrix Size:');
disp(size(behavioral_similarity));  % Expected: (n, n)

% Initialize correlation results per ROI
isc_behavior_corr = zeros(r, 1);

% Compute correlation for each ROI
for roi = 1:r
    % Extract ISC values for this ROI
    isc_values = squeeze(pairwise_temporal_ISC(:, :, roi));  % (n, n)

    % Vectorize both matrices (excluding diagonal to avoid self-correlation)
    isc_vector = isc_values(triu(true(n, n), 1));
    behavior_vector = behavioral_similarity(triu(true(n, n), 1));

    % Compute Pearson correlation between ISC and behavioral similarity
    isc_behavior_corr(roi) = corr(isc_vector, behavior_vector, 'rows', 'complete'); %Use rows complete to handle NaNs
end

% Display correlation results
disp('Correlation Between Pairwise ISC and Behavioral Similarity (per ROI):');
disp(isc_behavior_corr);

% Save the results
save(fullfile(data_folder,'representational_similarity_results.mat'), 'behavioral_similarity', 'isc_behavior_corr');

disp('Representational similarity analysis completed and results saved.');


% Results:
% The correlation between ISC and behavioral similarity for each ROI
% was computed and the results are as follows:
% - ROI 1: 0.0419
% - ROI 2: -0.1088
% - ROI 3: -0.0330
% - ROI 4: 0.0494
% - ROI 5: -0.1036
% 
% Interpretation:
% The correlation values are relatively low and in some cases negative, 
% suggesting little to no strong relationship between behavioral similarity 
% and ISC in these ROIs. This may imply that behavioral similarity does 
% not strongly influence the degree of neural synchronization in the regions 
% examined, or that other factors might be contributing to ISC in these 
% brain areas. Further analysis or additional data may be needed to 
% explore these findings in more depth.
