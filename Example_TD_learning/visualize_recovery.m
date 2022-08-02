function visualize_recovery(PRF, voxels, params, ncolumns)

arguments
    PRF
    voxels
    params
    ncolumns = 4
end

% calculate number of rows:
nparams = length(params);
nrows = nparams / ncolumns;

if mod(nparams, ncolumns) > 0
    nrows = nrows + 1;
end

% Params to P:
for ii = 1 : nparams
    P(ii).mu_alpha = params(ii, 1);
    P(ii).mu_eta = params(ii, 2);
end

samples = 200;

grid = PRF.U(1).grid;

d1 = grid.alpha;
d2 = grid.eta;
x_offset = -d1(1);
y_offset = -d2(1);
x_range= d1(2)-d1(1);
y_range= d2(2)-d2(1);

figure;

for vidx = 1 : length(voxels)

    a = subplot(nrows, ncolumns, vidx);
    z = cpm_draw_voxel(PRF, voxels(vidx), 'alpha', 'eta', '', samples, false);
    hold on;

    imagesc(z)
    contour(z)
    circles(P(vidx), x_offset, y_offset, x_range, y_range, samples)
    xlabel('alpha')
    ylabel('eta')
    
    x_labs = linspace(d1(1), d1(2), 7);
    x_ticks = (x_labs + x_offset) * samples / x_range;
    xticks(x_ticks)
    y_labs = linspace(d2(1), d2(2), 7);
    y_ticks = (y_labs + y_offset) * samples / y_range;
    yticks(y_ticks)
        
    xlim([x_ticks(1), x_ticks(end)]);
    ylim([y_ticks(1), y_ticks(end)]);

    yticklabels(num2cell(round(y_labs,2)))
    xticklabels(num2cell(round(x_labs, 2)))
    xtickangle(45)


    title(sprintf('CPM for\n alpha = %4.2f, eta = %4.2f', P(vidx).mu_alpha, P(vidx).mu_eta))
    set(a, 'dataAspectRatio', [1, 1, 1]);
end


end

function circles(P, x_offset, y_offset, x_range, y_range, samples)

   x = (P.mu_alpha + x_offset) * samples / x_range;
   y = (P.mu_eta + y_offset) * samples / y_range;
    viscircles([x y ], 1.25, 'Color', 'r', 'LineWidth', 1, 'EnhanceVisibility', true);

end
