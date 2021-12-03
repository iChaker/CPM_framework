function PRF = cpm_specify(SPM,options,y,XYZmm,U,receptive_field,observation_function,outpath,varargin)
% must give in options: name, TE  + other optional args in spm_prf_analyse

%SPM
SPM.swd = outpath;

% VOI
xY.y = y;
xY.XYZmm=XYZmm;
VOI.Y=y; % no eigen summary method
VOI.xY=xY;

p_names = fieldnames(U(1).grid);
% priors -------------------------------------------------

if isempty(receptive_field)
    receptive_field='cpm_RF_Gaussian';
end
 [~,~,~,pE,pC] = feval(receptive_field,p_names);


%options
if isempty(observation_function)
    observation_function='cpm_obv_identity';
end

[~,pE_obv,pC_obv] = feval(observation_function);
o_names= fieldnames(pE_obv);
for i=1:length(o_names)
    pE.(o_names{i}) = pE_obv.(o_names{i});
    pC.(o_names{i}) = pC_obv.(o_names{i});
end



options.pE = pE;
options.pC = pC;
options.model = 'spm_prf_response';
options.voxel_wise = true;
cpm.obv_fun=observation_function;
cpm.RF_fun=receptive_field;
try cpm.other=varargin; catch end
options.cpm=cpm;


PRF  = spm_prf_analyse('specify',SPM,VOI,U,options);

end

