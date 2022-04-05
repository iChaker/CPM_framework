function [trial_n, timestep_n, stimuli_n, S, stimuli, C, onsets] = cpm_prepare_experiment_data(filename, trial_n, fileOptions)
%PREPARE_EXPERIMENT_DATA Loads and prepares an event.tsv file in BIDS
%format.
% Requires matlabs table function
arguments
   filename char
   trial_n double = NaN
   fileOptions.response_cue string = 'Response_window_onset'
   fileOptions.response_stimulus string = 'Stimulus_onset'
   fileOptions.response_wealth string = 'New_wealth_onset'
   fileOptions.old_wealth string = 'old_wealth'
   fileOptions.new_wealth string = 'new_wealth'
   fileOptions.stimulus_id string = 'gFactIncr'
end

% Read in BIDS table in .tsv, specifying NaNs.

input_file = readtable(filename, 'FileType', 'text', 'TreatAsEmpty', 'n/a');

% Extract trials of interest
response_window = input_file(input_file.event_type == fileOptions.response_cue, :);
stimulus_onset = input_file(input_file.event_type == fileOptions.response_stimulus, :);
new_wealth = input_file(input_file.event_type == fileOptions.response_wealth, :);

% Assertion, that extract files are of same length!
assert (height(new_wealth) == height(response_window) && height(new_wealth) == height(stimulus_onset)); 

% For ease:
T0s = response_window.onset;
old_wealth = response_window.(fileOptions.old_wealth);
TSs = stimulus_onset.onset;
Ss = stimulus_onset.(fileOptions.stimulus_id);
TRs = new_wealth.onset;
new_wealth = new_wealth.(fileOptions.new_wealth);

timestep_n = 5;

if isnan(trial_n)
    trial_n = length(T0s);
end

stimuli_n = length(unique(Ss)) + 1;
stimuli = unique(Ss);

fprintf('number of trials: %i, number of stimuli %i \n', trial_n, stimuli_n);

onsets = [T0s(1: trial_n), TSs(1: trial_n), TRs(1: trial_n)]';
onsets = onsets(:);


% Constructing wealth matrix
C = zeros(trial_n, timestep_n);
C(:, [1:3]) = repmat(old_wealth(1: trial_n), 1, 3);
C(:, [4:5]) = repmat(new_wealth(1: trial_n), 1, 2);

S = zeros(trial_n, timestep_n, stimuli_n);

% Constructing stimulus identies:
for ii = 1 : trial_n
   S(ii, 3, find(stimuli == Ss(ii))) = 1;
   S(ii, 2, end) = 1;
end

end
