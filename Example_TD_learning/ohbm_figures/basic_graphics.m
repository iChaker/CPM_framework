addpath(genpath('/mnt/projects/CPM/'))

experimental_data = readtable('../tsvfile1.tsv', 'FileType', 'text', 'TreatAsEmpty', 'n/a');

%% Trial trajectory

wealth_idx = ~isnan(experimental_data.old_wealth);
filtered_data = experimental_data(wealth_idx, :);


fig1 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 600, 500]);

plot(filtered_data.onset, filtered_data.old_wealth, ':o', 'MarkerFaceColor', 'blue', 'MarkerSize', 6, 'LineWidth', 2)
% set(gca,'XColor','none', 'YColor', 'none')

xlabel('time')
ylabel('wealth')
fontname(fig1, 'Arial')
fontsize(fig1, scale=1.5)

exportgraphics(fig1,'wealth_trajectory.png','Resolution',600)
%% 
colormap parula

PRF = load('../simulationfiles/twotau_PRFn_estimated.mat');
PRF = PRF.PRFn{2};

fig2 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 1200, 1200]);

response = spm_prf_response(PRF.M.pE{1}, PRF.M, PRF.U, 'get_response', PRF.U(1).gridpoints);

response = reshape(response, 41, 41);

pcolor(response)
pbaspect([1 1 1])
set(gca,'XColor','none', 'YColor', 'none')
exportgraphics(fig2,'field.png','Resolution',600)

%%

bold = PRF.xY.y(:, 1);

fig3 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 600, 500]);

plot([1 : length(bold)] * 0.592, bold, '-', 'MarkerFaceColor', 'blue', 'LineWidth', 2)
% set(gca,'XColor','none', 'YColor', 'none')

xlabel('time')
ylabel('BOLD signal')
fontname(fig3, 'Arial')
fontsize(fig3, scale=1.5)

exportgraphics(fig3,'bold_trajectory.png','Resolution',600)

%%
[x, z ] =spm_prf_response(PRF.M.pE{1}, PRF.M, PRF.U);

z_u = z.u;
t = [1 : length(z_u)] * z.dt;

fig4 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 600, 500]);

plot(t, z_u, '-', 'MarkerFaceColor', 'blue', 'LineWidth', 2)
% set(gca,'XColor','none', 'YColor', 'none')

xlabel('time')
ylabel('neuronal signal')
fontname(fig4, 'Arial')
fontsize(fig4, scale=1.5)

exportgraphics(fig4,'u_trajectory.png','Resolution',600)

%%
colormap turbo

for ii = [1 : 50 : 504 ]
fig5 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 1200, 1200]);


response =PRF.U(ii).signals;

pcolor(response)
pbaspect([1 1 1])
set(gca,'XColor','none', 'YColor', 'none')
exportgraphics(fig5,['erp_response_' num2str(ii) '.png'],'Resolution',600)
end