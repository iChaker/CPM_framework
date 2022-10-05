function [PRFname, PRFn] = simulate_data(cfg_file, REDO, use_par)

arguments
    cfg_file = 'ideal_grid_cfg.json'
    REDO = true
    use_par =true
end


%%  ===========================================
% =================== SETUP ==================== =
% =============================================
% Random seed for reproducibility of noise etc.
rng(2022)
%% Creating directories to store simulation files 
cfg = spm_jsonread(cfg_file);

tmpdir = cfg.tmpdir;
simulationdir = cfg.simulationdir;
simulationname = cfg.simulationname;

mkdir(tmpdir); % Create a directory for temporary files which will be ignored.
mkdir(simulationdir);
%% ======== SIMULATION SETUP =======================


% We can freely chose for which TR we want to simulate data, we use the TR from
% a previous project of ours. 
TR = 0.592;
% The file containing the experimental inputs.
file="tsvfile1.tsv";
% dt is the micro resolution, which is often 8 or 20, but for example when slice-time
% correction has been produced should be set to the number of slices.
dt = TR / 28; % Assuming classical SPM approach 
% As this is a simulation, nscans are inferred from the data, but you can change
% it here manually:
nscans = nan;
% I see two possibilities in generating the data for our recovery study
% 1. Using a grid and moving the parameters to certain areas within the grid,
% where the grid is similar to the one we use to recover the values.
% 2. Using the point estimates (from U) and generating the BOLD signal directly
% from the RPE trajectory.


nnoise = length(cfg.sim_noise);
nnames = length(cfg.param_names);

p_grid = {};
rec_grid = {};
% These are hte options used by BayesPRF:
options = {};
options.name='simulation'; % PRF file tag
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
options.model='spm_prf_response'; % We use a modified version of this function.
%
rec_options = options;
rec_options.name = 'recovery';

mu_simulations = {};
sigma_simulations ={};


for pn = cfg.param_names(:)'
    p_grid.(pn{1}) = [cfg.(['p_grid_' pn{1}])(:)', cfg.p_resolution];
    rec_grid.(pn{1}) = [cfg.(['rec_grid_' pn{1}])(:)', cfg.rec_resolution];

    try 
        options.mu.(pn{1}) = cfg.(['p_mu_' pn{1}])(:)';
    catch
        disp("Mu boundaries not defined for simulation.")
    end
    
    try
        rec_options.mu.(pn{1}) = cfg.(['rec_mu_' pn{1}])(:)';
    catch
        disp("Mu boundaries not defined for recovery.")
    end

    mu_simulations{end + 1, 1} = cfg.(['sim_mu_' pn{1}])(:);
    sigma_simulations{end + 1, 1} = cfg.(['sim_sigma_' pn{1}])(:);
end

% The resolution of the different grids, generation grid will only be used when
mu_grid = cell(nnames, 1);
[mu_grid{1 : nnames}] = ndgrid(mu_simulations{:});

sigma_grid  = cell(nnames, 1);
[sigma_grid{1:  nnames}] = ndgrid(sigma_simulations{:});

for np = 1 : nnames
    mu_grid{np} = mu_grid{np}(:);
    sigma_grid{np} = sigma_grid{np}(:); 
end

nparams = length(mu_grid{1});
nstandard = length(sigma_grid{1});

%% =============== SIMULATING DATA ===============
% Preparing behavioural data for the simulation:
% Implementation note here: cpm_precompute requires data to only have the fields
% dur, dt, and ons. Everything else is further passed to the generation
% functions. In this case S and C which contain the complete serial compound
% vector and the wealth trajectory.
data = [];
[trial_n, ~, ~, S, ~, C, onsets] = cpm_prepare_experiment_data(file);
data.dur = ones(3 * trial_n, 1) .* 2 .* dt;
data.dt = ones(3 * trial_n, 1) .* dt;
data.ons = onsets;
data.C = C;
data.S = S;

% The BayesPRF requires an SPM struct, but only a few fields from there, which
% we generate here:
SPM = {};
SPM.xY.RT = TR;
SPM.swd = ''; % this field will be overwritten by CPM.

% As we are using a simulation we can infer the number of scans from the onsets
% and the TR:
if isnan(nscans)
    nscans = ceil(max(onsets) / TR);
end

%% Select grid model & precompute
model = "cpm_grid_RPE"; % user defined model, see cpm_grid_template for details

fixedparams.gamma=0.97;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.Lambda=1.0;

% The number of voxels is hte number of combinations and the number of noise
% levels:
nvoxels = nparams * nstandard * nnoise;
% Pre allocating:
y = zeros(nscans, nvoxels);
states = [];
% The voxel location is meaningless in our simulation, 
% but we used it to store noise level and indexes in the parameter space
xyz = zeros(3, nvoxels); 
xyz_def = cell(3, nvoxels);

inverse_fun = @(x, xmin, xmax) (norminv((x - xmin) ./ (xmax - xmin)));
%% 
if ~ exist(fullfile(simulationdir, [simulationname, '_simVOI.mat']), 'file') || REDO

    cc = 1;
    for sidx = 1 :  nstandard
        for vidx = 1 :  nparams
            % For safety clearing previously set values
            grid = [];
            U_voi = [];
            % BOLD generation block:
            % in the grid generation we use a grid over a number of points and
            % generat the BOLD response over the grid using the generative
            % process described by the CPM model
            if REDO && cc == 1 % the Grid needs to be oncly created once, so here we load it or recreat it if necessary
                U_voi = cpm_precompute(model, p_grid, fixedparams, data, fullfile(tmpdir, [simulationname '_simU_grid.mat']), true);
            else
                 U_voi = cpm_precompute(model, p_grid, fixedparams, data, fullfile(tmpdir, [simulationname '_simU_grid.mat']), false);
            end
            % We generate a temporary PRF file, to generate prior values and some additional options:
           tmpPRF = cpm_specify(SPM, options, zeros(nscans, 1), ...
                                                  zeros(3, 1), U_voi, 'cpm_RF_Gaussian', 'cpm_obv_int', [tmpdir filesep]);
            % We get the prior values of the latent(!!) parameters
            tmpPE = tmpPRF.M.pE{1};
            % As we define the to be simulated values in native space we have to
            % transform them into the laten space:
            
            pc = 1;
            for pn = cfg.param_names(:)'
                tmpPE.(['lmu_' pn{1}]) = inverse_fun(mu_grid{pc}(vidx),  tmpPRF.options.cpm.mu.(pn{1})(1), ...
                                                                                               tmpPRF.options.cpm.mu.(pn{1})(2));
                tmpPE.(['lsigma_' pn{1}])  =sigma_grid{pc}(sidx); 
                pc = pc + 1;
            end
            % we overwrite values accordingly in the prior structure:
            disp(cpm_get_true_parameters(tmpPE, tmpPRF.M, tmpPRF.U))
            [tmpy, stat] =  spm_prf_response(tmpPE, tmpPRF.M, tmpPRF.U);
            % Saving hidden states of the model for later investigation
            states = [states, stat.u];

            % After having generated the BOLD signal we add simulated noise:
            for nidx = 1 : nnoise
                % normrand with 0 variance returns the distributions mean:
                y(:, cc) = tmpy +  normrnd(0, cfg.sim_noise(nidx), [nscans,1]);
                xyz(:, cc) = [vidx; sidx; nidx];
                
                mus = zeros(1, nnames);
                sds = zeros(1, nnames);       
                for xx = 1 : nnames
                    mus(xx) = mu_grid{xx}(vidx);
                    sds(xx) = sigma_grid{xx}(sidx);
                end
         
                xyz_def(:, cc) = {mus; sds; cfg.sim_noise(nidx)};
                cc = cc + 1;
            end
        end
    
    end
    % Storing and saving the VOI in SPM format:
    VOI.xY.y = y;
    VOI.xY.XYZmm = xyz;
    VOI.xyz_def = xyz_def;
    % Mean values for completeness
    VOI.y = mean(y, 2); % for completeness
    save(fullfile(simulationdir, [simulationname, '_simVOI.mat']), 'VOI')

else
    VOI = load(fullfile(simulationdir, [simulationname '_simVOI.mat']));
    VOI = VOI.VOI;
end

%% ===================== MODEL INVERSION ======================
outpath=simulationdir;  % PRF file path
% Now we generate / precompute the recovery grid:

U_recovery = cpm_precompute(model, rec_grid, fixedparams, data, fullfile(simulationdir, [simulationname '_U_recovery']), REDO);

% And specify the PRF for recovery:
PRF = cpm_specify(SPM, options, VOI.xY.y, VOI.xY.XYZmm, ...
                   U_recovery, 'cpm_RF_Gaussian', 'cpm_obv_int', outpath);

voxels = []; % or select subset of voxels of interest
PRF.M.noprint = 0; % to suppress spm outputs
PRFname = fullfile(simulationdir, [simulationname '_PRFn.mat']);


if ~exist(PRFname, 'file') || REDO
    PRFn = cpm_estimate(PRF, voxels, use_par);
    save(PRFname, 'PRFn'),
else
    PRFn = load(PRFname);
    PRFn = PRFn.PRFn;
end

end