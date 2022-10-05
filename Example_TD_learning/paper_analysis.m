% Clear console etc.
clear;
clc;
close all;
% Adding paths
addpath('utils_rpe')
addpath('utils_simulation')
addpath(genpath('../CPM_Toolbox'))
addpath(genpath('../BayespRF-master'))
% Adding toolboxes, check that you have the right path here!
addpath('/mnt/projects/LogUtil/Neuroimaging_RewardCoding/scripts/cpm.bu/toolboxes/spm12/')
BASEDIR = pwd;
cd('/mnt/projects/LogUtil/Neuroimaging_RewardCoding/scripts/cpm.bu/toolboxes/VBA-toolbox/')
VBA_setup
cd(BASEDIR);


[prf_path, PRF] = simulate_data('ideal_grid_cfg.json');


