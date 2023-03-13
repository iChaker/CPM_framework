clear;
close all;
clc;

rng(1,'philox')

addpath('/mnt/projects/CPM/toolboxes/spm12/')
addpath('utils_simulation')
addpath(genpath('../CPM_Toolbox'))
addpath(genpath('../BayespRF-master'))

%
M = {};
M.l = 1;
M.dt = 0.0211;
M.ns = 2690;
M.f = 'spm_prf_fx';
M.g = 'spm_prf_gx';
M.m = 1; 
M.x = [0, 0 ,0 , 0];

P = {};
P.transit =  0;
P.decay = 0;
P.epsilon =  -0.7800;
%
no_sims = 5;

z = zeros(1, M.ns  * 28);
non_spar = zeros(1, 28);


times = zeros(no_sims, size(non_spar, 2));
times_sparse = zeros(no_sims, size(non_spar, 2));

signalout = zeros(no_sims, size(non_spar ,2), M.ns);
signalout_sparse = zeros(no_sims, size(non_spar ,2), M.ns);

length_u = M.ns * 28 - 32;

for ii = 1 :  28
    
    tmp_z = z;

    for kk = 1 : ii
        nsteps = length(kk : 28 : length_u);
        tmp_z(kk : 28 : length_u) = rand(1, nsteps);
    end

    non_spar(ii) = length(tmp_z(tmp_z ~= 0)) / (M.ns  * 28);

    for jj = 1 : no_sims
        tic
        ty = spm_int_sparse(P, M, tmp_z');
        times_sparse(jj, ii) = toc;
        signalout_sparse(jj, ii, :) = ty;
            
        tic
        ty = spm_int(P, M, tmp_z');
        times(jj, ii) = toc;
        signalout(jj, ii, :) = ty;
    end
end

%%
fig = figure('Position', [0, 0, 1400, 400]); 


% limits:
y_lim = [0, max([max(times(:), max(times_sparse(:)))])];
y_lim(1)  = y_lim(1) + 0.25;
non_spar_label = {};

for ii = 1  : length(non_spar)
    if mod(ii, 2) == true
        non_spar_label{ii} = [num2str(round(non_spar(ii) * 100, 0)) '%']; 
    end
end

subplot(1, 3, 1); 
boxchart(times); 
xticklabels(non_spar_label)
ylabel('time in s')
xlabel('non-zero values (percent)')
ylim(y_lim)
title('spm\_int')

subplot(1, 3, 2); 
boxchart(times_sparse); 
xticklabels(non_spar_label)
ylabel('time in s')
xlabel('non-zero values (percent)')
title('spm\_int\_sparse')
ylim(y_lim)

subplot(1, 3, 3);
plot(round(non_spar * 100, 1), mean(times) ./ mean(times_sparse))
xticks(round(non_spar * 100, 1))
xticklabels(non_spar_label)
xlabel('non-zero values (percent)')
ylabel('Speed up factor')
title(['Speed up'])

sgtitle('Optimization of spm\_int')

cpm_savefig(fig, 'results/figs3_speed_up.png')