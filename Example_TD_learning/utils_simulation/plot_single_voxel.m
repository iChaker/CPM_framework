

Ep = PRF.Ep{1};

plot_names = {'eta', 'tau'};

transformers = {@cpm_logit_inv, []};
transform_names = {'alpha', []};

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
XYZ = [];

for ii = 1 : length(param_names)
XYZ = [XYZ, xyz{ii}(:)];
end

%%
[~, z] = cpm_prf_get_ppd(PRF, XYZ, 1, 100, true);

z = reshape(z, resolution{:});

plt_idx = nan(1, length(plot_names));

for cc = 1 : length(plot_names)
    plt_idx(cc) = find(contains(param_names, plot_names{cc}));
end

zdims = 1 : length(param_names);
zp = permute(z, [plt_idx, zdims(~ismember(zdims, plt_idx))]);

param_names = param_names( [plt_idx, zdims(~ismember(zdims, plt_idx))]);

true_ep = cpm_get_true_parameters(PRF, 1);
true_pe = cpm_get_true_parameters(PRF.M.pE{1}, PRF.M, PRF.U);
% A is given, up to 8 dims here.

cols = repmat({':'},1,ndims(zp));
for ix = length(plot_names) + 1 : length(param_names)
    % slice at posterior estimate of parameter
     [~, post] = min(abs(griddim.(param_names{ix}) - true_ep.(['mu_'  param_names{ix}])));
    cols{ix} = post;
end

zp_slc = zp(cols{:});

figure;

if length(plot_names) == 1
    plot(zp_slc)
else
    
        hold on
        rec_pos = [0 - diff(pmus.(plot_names{1})) / 2, 0 - diff(pmus.(plot_names{2})) / 2, ...
                            diff(pmus.(plot_names{1})) , diff(pmus.(plot_names{2}))];

        imagesc(griddim.(plot_names{1}), griddim.(plot_names{2}), zp_slc)

        contour(griddim.(plot_names{1}), griddim.(plot_names{2}), zp_slc)
       
        rectangle('Position', rec_pos, 'EdgeColor','white', 'LineWidth', 2, 'LineStyle', '-');

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
    
        x_post = true_ep.(['mu_' plot_names{1}]) + 2 *  true_ep.(['sigma_' plot_names{1}]) * cos(t);
        y_post = true_ep.(['mu_' plot_names{2}]) +  2 *  true_ep.(['sigma_' plot_names{2}]) * sin(t);
        
        plot(x_post, y_post, 'white', 'LineWidth', 2, 'LineStyle', '--')
        

    pbaspect([1 1 1])
    %daspect([1 1 1])
    set(gca, 'LooseInset', get(gca,'TightInset'))
end


