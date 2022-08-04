function fig = visualize_recovery(PRF, voxels, params, ncolumns, posterior, field_resolution, param_names, ppd_samples)

arguments
    PRF
    voxels
    params
    ncolumns = 4
    posterior = true
    field_resolution = 100
    param_names = {'alpha', 'eta'}
    ppd_samples = 500
end

% calculate number of rows:
nparams = length(params);
nrows = nparams / ncolumns;

if mod(nparams, ncolumns) > 0
    nrows = nrows + 1;
end

% Params to P:
for ii = 1 : nparams
    P(ii).(['mu_' param_names{1}]) = params(ii, 1);
    P(ii).(['mu_' param_names{2}]) = params(ii, 2);
end

grid = PRF.U(1).grid;

x_bins = linspace(grid.(param_names{1})(1), grid.(param_names{1})(2), field_resolution);
y_bins = linspace(grid.(param_names{2})(1), grid.(param_names{2})(2), field_resolution);
[x2,y2] = meshgrid(x_bins,y_bins);
xy = [x2(:) y2(:)];

d1 = grid.(param_names{1}); 
d2 = grid.(param_names{2});
x_offset = -d1(1);
y_offset = -d2(1);
x_range= d1(2)-d1(1);
y_range= d2(2)-d2(1);

fig = figure('Color', 'none', 'Units', 'pixels', 'Position', [0, 0, 1600, 1600]);

for vidx = 1 : length(voxels)

    a = subplot(nrows, ncolumns, vidx);
    [~, ~, zprior, zpost] = spm_prf_get_ppd(PRF, xy, voxels(vidx), ppd_samples);
    
    if posterior
        z = zpost;
    else
        z = zprior;
    end
    hold on;
    z = reshape(z, field_resolution, field_resolution);
    imagesc(z)
    contour(z)
    circles(P, x_offset, y_offset, x_range, y_range, field_resolution, param_names)
    xlabel('alpha')
    ylabel('eta')
    
    x_labs = linspace(d1(1), d1(2), 7);
    x_ticks = (x_labs + x_offset) * field_resolution / x_range;
    xticks(x_ticks)
    y_labs = linspace(d2(1), d2(2), 7);
    y_ticks = (y_labs + y_offset) * field_resolution / y_range;
    yticks(y_ticks)
        


    yticklabels(num2cell(round(y_labs,2)))
    xticklabels(num2cell(round(x_labs, 2)))
    xtickangle(45)
    
    xlim([x_ticks(1), x_ticks(end)]);
    ylim([y_ticks(1), y_ticks(end)]);

    title(sprintf('CPM for\n alpha = %4.2f, eta = %4.2f', P(vidx).mu_alpha, P(vidx).mu_eta))
    set(a, 'dataAspectRatio', [1, 1, 1]);
end


end

function circles(P, x_offset, y_offset, x_range, y_range, samples, param_names)
   
for ii = 1 : length(P)
       x = (P(ii).(['mu_' param_names{1}]) + x_offset) * samples / x_range;
       y = (P(ii).(['mu_' param_names{2}]) + y_offset) * samples / y_range;
        plot(x, y, 'o', 'Color', 'white')
end
end
