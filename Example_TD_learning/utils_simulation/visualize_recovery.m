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

height = 1200;
fig = figure('Color', 'none', 'Units', 'pixels', 'Position', [1, 1, 3/3 * height, height]);

tiles = tiledlayout(nrows, ncolumns); 

for vidx = 1 : length(voxels)
    nexttile
    %a = subplot(nrows, ncolumns, vidx);
    [~, z] = cpm_prf_get_ppd(PRF, xy, voxels(vidx), ppd_samples, posterior);
    
    hold on;
    z = reshape(z, field_resolution, field_resolution);
    imagesc(z)
    contour(z)
    circles(P, x_offset, y_offset, x_range, y_range, field_resolution, param_names, vidx)

    
    x_labs = linspace(d1(1), d1(2), 7);
    x_labs(1) = x_labs(1) + 0.001;
    
    x_ticks = (x_labs + x_offset) * field_resolution / x_range;
    x_ticks(1) = 0.5;
    xticks(x_ticks)
    y_labs = linspace(d2(1), d2(2), 7);
    y_ticks = (y_labs + y_offset) * field_resolution / y_range;
    y_ticks(1) = 0.5; 
    yticks(y_ticks)

    
    yticklabels(num2cell(round(y_labs,2)))
    xticklabels(num2cell(round(x_labs, 2)))
    xtickangle(45)
    
    a = gca;
%     set(a, 'dataAspectRatio', [1, 1, 1]);
    ttl = sprintf('CPM for\n%s = %4.2f\n%s = %4.2f', param_names{1}, P(vidx).(['mu_' param_names{1}]),  ...
                                                                                  param_names{2}, P(vidx).(['mu_' param_names{2}]));
    %title(ttl);
    pbaspect([1 1 1])
    daspect([1 1 1])
    set(gca, 'LooseInset', get(gca,'TightInset'))
    
    xlim([x_ticks(1), x_ticks(end)]);
    ylim([y_ticks(1), y_ticks(end)]);
end
tiles.TileSpacing = 'none';
tiles.Padding = 'none';

    xlabel(tiles, param_names{1})
    ylabel(tiles, param_names{2})
    %tiles.Position(3)=.999;
    
    
end

function circles(P, x_offset, y_offset, x_range, y_range, samples, param_names, idx)
   ax = gca;
   mcolor ='white'; % [174,1,126] ./ 255;
for ii = 1 : length(P)
       x = (P(ii).(['mu_' param_names{1}]) + x_offset) * samples / x_range;

       if x > max(ax.XLim)
           x = max(ax.XLim) - 0.5;
       elseif x < min(ax.XLim)
           x = min(ax.XLim) + 0.5;
       end
       

       y = (P(ii).(['mu_' param_names{2}]) + y_offset) * samples / y_range;       
       
        if y > max(ax.YLim)
           y = max(ax.YLim) -1;
       elseif y < min(ax.YLim)
           y = min(ax.YLim) +1;
        end

        if ii == idx
            plot(x, y, '+', 'Color', mcolor, 'MarkerFaceColor',  mcolor,  'MarkerSize', 8, 'LineWidth', 2)
        else
            plot(x, y, 'o', 'Color', mcolor , 'MarkerSize', 4, 'LineWidth', 1)
        end

end
end
