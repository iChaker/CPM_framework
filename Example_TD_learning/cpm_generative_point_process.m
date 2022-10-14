function [y] = cpm_generative_point_process(params, model, fixedparams, data, TR, nscans, tmpdir, simulationname, observation_function)
%DATA_GENERATIVE_PROCESS Summary of this function goes here
% Generate data for a set of parameters in a given grid. 
% Full grid has to be provided and can be reduced  / collapsed at a given point for using CPM. 

arguments
    params
    model
    fixedparams
    data
    TR
    nscans
    tmpdir = ''
    simulationname = 'point'
    observation_function = 'cpm_obv_int'
end


pnames = fieldnames(params);
grid = {};

for pn = pnames'
    grid.(pn{1}) = [params.(pn{1}), params.(pn{1}), 1];
end

% Generate a point grid:
U = cpm_precompute(model, grid, fixedparams, data, fullfile(tmpdir, [simulationname '_simU_point.mat']), true);
[y, ~, ~] = cpm_generate_from_U(U, 1, TR, nscans, observation_function);


end

