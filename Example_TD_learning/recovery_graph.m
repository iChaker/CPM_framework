
% load sims and PRFs by hand
b = 100;

grid = PRF.U(1).grid;
d1 = [0 1 grid.tau(3)];
d2 = grid.eta;
x_offset = 0;
y_offset = -d2(1);
x_range= d1(2)-d1(1);
y_range= d2(2)-d2(1);
params = cat(1,sim_nul.xY.Params{:},sim_taueta.xY.Params{:},sim_tau.xY.Params{:});


figure;
t = tiledlayout(3,5);
xlabel(t,'tau');
ylabel(t,'eta');
sgtitle("Posterior distributions");

draw_alpha(PRF_nul,1,'tau','eta','',b,6);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);

draw_alpha(PRF_taueta,1,'tau','eta','',b,2);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
draw_alpha(PRF_taueta,2,'tau','eta','',b,3);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
draw_alpha(PRF_taueta,3,'tau','eta','',b,4);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
draw_alpha(PRF_taueta,4,'tau','eta','',b,5);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);


cpm_draw_voxel(PRF_tau,1,'tau','eta','',b,7);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
cpm_draw_voxel(PRF_tau,2,'tau','eta','',b,8);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
cpm_draw_voxel(PRF_tau,3,'tau','eta','',b,9);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
cpm_draw_voxel(PRF_tau,4,'tau','eta','',b,10);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);

cpm_draw_voxel(PRF_taueta,5,'tau','eta','',b,12);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
cpm_draw_voxel(PRF_taueta,6,'tau','eta','',b,13);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
cpm_draw_voxel(PRF_taueta,7,'tau','eta','',b,14);
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);
cpm_draw_voxel(PRF_taueta,8,'tau','eta','',b,15 );
circles(params,PRF,x_offset,y_offset,x_range,y_range,b);





function circles(params,PRF,x_offset,y_offset,x_range,y_range,b)

    for i=1:length(params)
       P = params(i);
       P = cpm_get_true_parameters(P,PRF.M,PRF.U);
       alpha = cpm_logit_inv(P.mu_tau);
       x = (alpha+x_offset) * b / x_range;
       y = (P.mu_eta+y_offset) * b / y_range;

        viscircles([x y ], 1.25,'Color','r','LineWidth',1,'EnhanceVisibility', false);
    end

end




