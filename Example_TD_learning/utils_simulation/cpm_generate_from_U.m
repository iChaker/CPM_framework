function [y, z, grid_z] = cpm_generate_from_U(U, idx, TR, nscans, observation_function)
% FUNCTION NAME
% cpm_generate_from_U
% 
% DESCRIPTION
% Generate data using the observation functions of the CPM/PRF toolbox,
% using only the raw inputs (U) from the model, interpreting them directly
% as the input "z". An idx has to be provided
% to indicate the space in the grid.
% It also allows the use of different TRs and could serve in model
% development. 
%
% INPUT:
%   U - the trial information
%   idx - idx for the set of parameters in U(:).signals1D
%   TR - the TR for which the bold signal should be generated.
%   nscans - number of datapoints to be generated, optional
%   observation_function - observation functon of the model, optional,
%   defaults to cpm_obv_int.
% 
% OUTPUT:
%   y - generated observed (i.e. BOLD) signal 
%   z - transformed input from U into continuous micro states.
%   grid_z - Signals at idx, from U.
%
% ASSUMPTIONS AND LIMITATIONS:
% Currently, is very rigid and only for use with cpm_obv_int.
%
% REVISION HISTORY:
% 29/06/2022 - SRSteinkamp
%   * First implementation. 

if nargin < 4
    nscans = round(max([U(:).ons]) ./ TR) + 1;
end

if nargin < 3
    TR = 0.592;
end

if nargin < 5
    observation_function = str2func('cpm_obv_int');
else
    observation_function = str2func(observation_function);
end

M = struct();

M.x = [0;0;0;0];
M.f = 'spm_prf_fx';
M.g = 'spm_prf_gx';
M.n = 4;
M.u = 1;
M.l = 1; 
M.m = 1;
M.u = 0;
M.ns = nscans;
M.dt = U(1).dt;

grid_z = [U(:).signals1D];

grid_z = grid_z(idx, :);

[~, pE] = cpm_obv_int();

% Adding info to U
bins_per_second = 1 / U(1).dt;
bins_per_TR     = bins_per_second * TR;

if ~floor(bins_per_TR)==bins_per_TR
    
    warning('Bins per TR is not an integer value, please check dt and TR input values');
    bins_per_TR = round(bins_per_TR);

elseif mod(bins_per_TR, 1) ~= 0

    warning('Bins per TR is not an integer value, please check dt and TR input values');
    bins_per_TR = round(bins_per_TR);
    
end


for t = 1:length(U)
    start_bin = ceil( U(t).ons / U(t).dt) + 1;
    end_bin   = start_bin + (U(t).dur * bins_per_second) - 1;

    U(t).ind   = start_bin : end_bin;           % Microtime bins (from 1)
    U(t).nbins = nscans * bins_per_TR;            % Total bins
end     

% From u to z:


 % Number of volumes = inputs
n  = length(U); 
nbins = max(U(1).nbins, max(U(end).ind));

z = zeros(1,nbins);
for t = 1 : n    
    % Microtime index for this volume
    ind = U(t).ind;  

    resp = U(t).signals1D(idx)';
    z(ind) = z(ind) + sum(resp);

end

Z.u = z';
Z.dt = M.dt;

y = observation_function(pE, M, Z);

end