function [trial_n, timestep_n, stimuli_n, S, stimuli, C, onsets] = cpm_prepare_experiment_data(filename)
%PREPARE_EXPERIMENT_DATA Loads and prepares an event.tsv file in BIDS
%format.
% Requires matlabs table function

% Read in BIDS table in .tsv, specifying NaNs.
input_file = readtable(filename, 'FileType', 'text', 'TreatAsEmpty', 'n/a');

% Extract trials of interest
response_window = input_file(input_file.event_type == "Response_window_onset", :);
stimulus_onset = input_file(input_file.event_type == "Stimulus_onset", :);
new_wealth = input_file(input_file.event_type == "New_wealth_onset", :);

% Assertion, that extract files are of same length!
assert (height(new_wealth) == height(response_window) && height(new_wealth) == height(stimulus_onset)); 
 
% For ease:
T0s = response_window.onset;
old_wealth = response_window.old_wealth;
TSs = stimulus_onset.onset;
Ss = stimulus_onset.gFactIncr;
TRs = new_wealth.onset;
new_wealth = new_wealth.new_wealth;
Rs = new_wealth - old_wealth;

timestep_n = 5;
trial_n = length(T0s);

stimuli_n = length(unique(Ss)) + 1;
stimuli = unique(Ss);

fprintf('number of trials: %i, number of stimuli %i \n', trial_n, stimuli_n);

onsets = [T0s, TSs, TRs]';
onsets = onsets(:);

S = zeros(trial_n, timestep_n, stimuli_n);

% Constructing wealth matrix
C = zeros(trial_n, timestep_n);
C(:, [1:3]) = repmat(old_wealth, 1, 3);
C(:, [4:5]) = repmat(new_wealth, 1, 2);
% Constructing reward
R = zeros(trial_n, timestep_n);
R(:, 4) = Rs;
stimuli_index = {};

% Constructing stimulus identies:
for ii = 1 : trial_n
   S(ii, 3, find(stimuli == Ss(ii))) = 1;
   S(ii, 2, end) = 1;
end

end
