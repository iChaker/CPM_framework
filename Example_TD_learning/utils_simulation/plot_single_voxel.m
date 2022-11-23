function plot_single_voxel(PRF, idx, plot_names, transformers, transform_names, num_samples, posterior)

arguments
    PRF
    idx = 1
    plot_names = {'tau', 'eta'}
    transformers = {@cpm_logit_inv, []}
    transform_names = {'alpha', 'eta'}
    num_samples = 100
    posterior = true
end


Ep = PRF.Ep{idx};

for ii = 1 : length(transformers)
    if isempty(transformers{ii})
        transformers{ii} = @(x) x;
        transform_names{ii} = plot_names{ii};
    end
        
end

if length(plot_names) > 2
    error("Not implemented!")
end

% Recover grid:
grid = {};
psigs = {};
pmus = {};
griddim = {};
griddim_noname = {};
resolution = {};

cc = 1;

param_names = fieldnames(PRF.U(1).grid);

for pn = param_names(:)'
    grid.(pn{1}) = PRF.U(1).grid.(pn{1});
    psigs.(pn{1}) = PRF.options.cpm.sigma.(pn{1});
    pmus.(pn{1}) = PRF.options.cpm.mu.(pn{1});
    griddim.(pn{1}) = linspace(grid.(pn{1})(1), grid.(pn{1})(2), grid.(pn{1})(3));
    griddim_noname{cc} = griddim.(pn{1});
    resolution{cc} = grid.(pn{1})(3);
    cc = cc + 1;
end

xyz = {};
[xyz{1 : length(param_names)}] = meshgrid(griddim_noname{:});

if length(griddim_noname) < 2
   xyz{1} =  xyz{1}(1, :);
end
    
XYZ = [];

for ii = 1 : length(param_names)
    XYZ = [XYZ, xyz{ii}(:)];
end

%%
[~, z] = cpm_prf_get_ppd(PRF, XYZ, idx, num_samples, posterior);

if length(resolution) < 2
    resolution{end + 1} = [];
end
z = reshape(z, resolution{:});

plt_idx = nan(1, length(plot_names));

for cc = 1 : length(plot_names)
    plt_idx(cc) = find(contains(param_names, plot_names{cc}));
end

zdims = 1 : length(param_names);

if length(zdims) > 1
    zp = permute(z, [plt_idx, zdims(~ismember(zdims, plt_idx))]);
else 
    zp = z;
end

param_names = param_names( [plt_idx, zdims(~ismember(zdims, plt_idx))]);
% transform_names = transform_names( [plt_idx, zdims(~ismember(zdims, plt_idx))]);
% transformers = transformers( [plt_idx, zdims(~ismember(zdims, plt_idx))]);

true_ep = cpm_get_true_parameters(PRF, idx);
% true_pe = cpm_get_true_parameters(PRF.M.pE{idx}, PRF.M, PRF.U);
% A is given, up to 8 dims here.

cols = repmat({':'},1,ndims(zp));
for ix = length(plot_names) + 1 : length(param_names)
    % slice at posterior estimate of parameter
     [~, post] = min(abs(griddim.(param_names{ix}) - true_ep.(['mu_'  param_names{ix}])));
    cols{ix} = post;
end

zp_slc = zp(cols{:});

if length(plot_names) == 1
    plot(zp_slc)
else
    
        hold on
        
        if length(zdims) < 2
            zp_slc = zp_slc * zp_slc';
        end
        pmus.(plot_names{1}) = transformers{1}(pmus.(plot_names{1}));
        pmus.(plot_names{2}) = transformers{2}(pmus.(plot_names{2}));
        
        rec_pos = [sum(pmus.(plot_names{1})) / 2 - diff(pmus.(plot_names{1})) / 2, ...
                          sum(pmus.(plot_names{2})) / 2 - diff(pmus.(plot_names{2})) / 2, ...
                            diff(pmus.(plot_names{1})) , ...
                            diff(pmus.(plot_names{2}))];
        
%         rec_pos(1 : 2 : end) = transformers{1}(rec_pos(1:2:end));
%         rec_pos(2 : 2 : end) = transformers{2}(rec_pos(2:2:end));
        
        pplot = pcolor(transformers{1}(griddim.(plot_names{1})), transformers{2}(griddim.(plot_names{2})), zp_slc);
        %contour(transformers{1}(griddim.(plot_names{1})), transformers{2}(griddim.(plot_names{2})), zp_slc)
        pplot.LineWidth = 0.5;
        pplot.FaceColor = 'interp';
        rectangle('Position', rec_pos, 'EdgeColor','white', 'LineWidth', 1, 'LineStyle', '--');

        t = -pi:0.01:pi;
%         x_min = true_pe.(['mu_' plot_names{1}]) +  2 * psigs.(plot_names{1})(1) * cos(t);
%         y_min = true_pe.(['mu_' plot_names{2}]) + 2 * psigs.(plot_names{2})(1) * sin(t);
%         x_max = true_pe.(['mu_' plot_names{1}]) + 2 *  psigs.(plot_names{1})(2) * cos(t);
%         y_max = true_pe.(['mu_' plot_names{2}]) + 2 * psigs.(plot_names{2})(2) * sin(t);
%         x_prior = true_pe.(['mu_' plot_names{1}]) +  2*  true_pe.(['sigma_' plot_names{1}]) * cos(t);
%         y_prior = true_pe.(['mu_' plot_names{2}]) +   2 * true_pe.(['sigma_' plot_names{2}]) * sin(t);
%         plot(x_min, y_min)
%         plot(x_max, y_max)
%         plot(x_prior, y_prior)
    
%         x_post = true_ep.(['mu_' plot_names{1}]) + 2 *  true_ep.(['sigma_' plot_names{1}]) * cos(t);
%         y_post = true_ep.(['mu_' plot_names{2}]) +  2 *  true_ep.(['sigma_' plot_names{2}]) * sin(t);
%         
%         x_post = transformers{1}(x_post);
%         y_post = transformers{2}(y_post);
       %plot(x_post, y_post, 'white', 'LineWidth', 2, 'LineStyle', '--')
        

    pbaspect([1 1 1])
    %daspect([1 1 1])
    set(gca, 'LooseInset', get(gca,'TightInset'))
end

end


