function visualize_recovery(PRF, voxels, params)

% Params to P:
for ii = 1 : length(params)
    P(ii).mu_alpha = params(ii, 1);
    P(ii).mu_eta = params(ii, 2);
end

samples = 200;

grid = PRF.U(1).grid;

d1 = grid.alpha;
% d1 = [0, 1, grid.tau(3)];
d2 = grid.eta;
x_offset = -d1(1);
y_offset = -d2(1);
x_range= d1(2)-d1(1);
y_range= d2(2)-d2(1);

figure;

for vidx = 1 : length(voxels)

    subplot(3, 4, vidx)
    cpm_draw_voxel(PRF, voxels(vidx), 'alpha', 'eta', '', samples, voxels(vidx));
    circles(P(vidx), PRF, x_offset, y_offset, x_range, y_range, samples)
end


end

function circles(params,PRF,x_offset,y_offset,x_range,y_range, samples)

    for i=1:length(params)
       P = params(i);
        %   P = cpm_get_true_parameters(P,PRF.M,PRF.U);
%        alpha = cpm_logit_inv(P.mu_tau);
        alpha = P.mu_alpha;
       x = (alpha+x_offset) * samples / x_range;
       y = (P.mu_eta+y_offset) * samples / y_range;
        viscircles([x y ], 1.25,'Color','r','LineWidth',2,'EnhanceVisibility', true);
    end

end
