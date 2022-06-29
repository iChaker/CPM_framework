function [true_P, tm_names, ts_names] = cpm_get_true_parameters(varargin)
% cpm_get_true_parameters 
%   This function transforms latent neuronal parameters. Returns the correct 
%   population field parameters for interpretation.
%
%   In: PRF and voxelnumber, alternatively, a structure containing all VL parameters, 
%       a PRF.M structure and a population response U 
%   
%   out: true population field parameters for that voxel

if nargin==2
    PRF=varargin{1};
    voxel=varargin{2};
    
    P = PRF.Ep{1,voxel};
    M = PRF.M;
    U = PRF.U;
    
elseif nargin==3
    P =varargin{1};
    M =varargin{2};
    U =varargin{3};

     
end
    

receptive_field = M.cpm.RF_fun;
[ ~, fn_true_moments,fn_true_scaling, fn_latent_moment_priors,fn_latent_scaling_priors ] = feval(receptive_field);
grid = U(1).grid;
p_names = fieldnames(grid);
[pE,~] = fn_latent_moment_priors();
m_names = fieldnames(pE);
[pE,~] = fn_latent_scaling_priors();
s_names = fieldnames(pE);

true_P= struct();

cc = 1;
L_names = {};

% Add latent parameters
for i=1:length(p_names)
   latent_moments = struct();
   for j=1:length(m_names)
       VL_latent = [ m_names{j} '_' p_names{i} ];
       latent_moments.(m_names{j})=P.(VL_latent);
       L_names{cc} = VL_latent;
       cc = cc + 1;
   end
   minp = U(1).grid.(p_names{i})(1);
   maxp = U(1).grid.(p_names{i})(2);
   true_moments = fn_true_moments(latent_moments,minp,maxp);
   tm_names = fieldnames( true_moments );
   for j=1:length(tm_names)
       VL_true = [ tm_names{j} '_' p_names{i} ];
       true_P.(VL_true) = true_moments.(tm_names{j});
   end
   
end

% Add scaling parameters
latent_scaling = struct();
for i=1:length(s_names)    
    latent_scaling.(s_names{i}) = P.(s_names{i});
end
true_scaling = fn_true_scaling(latent_scaling);
ts_names = fieldnames(true_scaling);
for i=1:length(ts_names)
   true_P.(ts_names{i}) =  true_scaling.(ts_names{i});
end

% Add non transformed parameters
all_names = fieldnames(P);

for i = 1 : length(all_names)
    if ~any(strcmpi(all_names{i}, [L_names, s_names]))
       true_P.(all_names{i}) = P.(all_names{i}); 
    end
end

% Sort true parameters to have same order as native
out_names = fieldnames(true_P);

sort_vec = zeros(length(out_names), 1);

% Works because, latent parameters (in all_names) have only a 'l' added
for i = 1 : length(out_names)
    sort_vec = sort_vec + (contains(all_names, out_names{i}) .* i);
end

true_P = orderfields(true_P, out_names(sort_vec));

end

