function U = cpm_precompute(model,grid,fixedparams,data,precomp_file)
% cpm_precompute Summary of this function goes here
%   INPUTS:
%       -model: name computational model of interest. The user would 
%               implement this function. See cpm_grid_template for details.
%       -grid: a structure that defines the grid over the parameter space
%       (grid.parami = [parami_min parami_max N_points]. The fieldnames
%       should correspond with the parameters used in the computational model function
%       -fixedparams: hyperparameters that are passed to the computational model
%       -data: stimulus data to be used for computing neuronal signals. See
%       cpm_grid_template for details
%       -precomp_file: filename of the resulting population field U
%
%   OUTPUTS:
%       -U: population response containing all possible responses.



if exist(precomp_file, 'file') == 2
    disp(['precomputation available: ' precomp_file])
    load(precomp_file);  

else

    linspaces = {};
    p_names = fieldnames(grid);

    for i=1:numel(p_names)
       dims(i)= grid.(p_names{i})(3);
       linspaces{i} = linspace(grid.(p_names{i})(1),...
                               grid.(p_names{i})(2),...
                               grid.(p_names{i})(3));
       idxs = linspace(1,...
                       grid.(p_names{i})(3),...
                       grid.(p_names{i})(3));
       idxcell{i}=idxs;   
                           
    end


    ons = data.ons;
    argmatrix = combvec(linspaces{:});
    argmatrix = argmatrix';

    idxmatrix = combvec(idxcell{:});
    idxmatrix = idxmatrix';


    p_values = num2cell(argmatrix(1,:),1);
    freeparams = cell2struct(p_values,p_names,2);
    tt = tic;
    signaltmp = feval(model,freeparams,fixedparams,data);
    tt=toc(tt);
    
    disp(' PRECOMPUTING MESHGRID USING GRID: ');
    disp(grid);
    disp(' ETA :');
    disp(duration(0,0,tt*size(idxmatrix,1)));
    
    precomputed = zeros([length(signaltmp),dims]);
    
    
    for i=1:size(idxmatrix,1)

        p_values = num2cell(argmatrix(i,:),1);
        freeparams = cell2struct(p_values,p_names,2);
        
        signal = feval(model,freeparams,fixedparams,data);

        idx = num2cell(idxmatrix(i,:));
        precomputed(:,idx{:})= signal; 
    end
    
    microtime = true;
    U(length(signal)) = struct();
    S= {};
    S.subs = repmat({':'},1,ndims(precomputed));
    S.type = '()';

    for nt = 1 : length(signal)
        S.subs{1}=nt;
        signals=squeeze(subsref(precomputed,S));
        
        signals1D = zeros(length(idxmatrix),1);
        for i=1:length(idxmatrix)
            idx = num2cell(idxmatrix(i,:));
            signals1D(i)=signals(idx{:});
        end
        
    
        U(nt).signals = signals;
        U(nt).signals1D = signals1D;
        try
            try 
                U(nt).ons = data.ons(nt);
                U(nt).dur = data.dur(nt);
                U(nt).dt = data.dt(nt);

            catch
                U(nt).ons = data(nt).ons;
                U(nt).dur = data(nt).dur;
                U(nt).dt = data(nt).dt;
            end
        catch
            try
                U(nt).ons = data.ons(1);
                U(nt).dur = data.dur(1);
                U(nt).dt = data.dt(1); 
                microtime  = false;
            catch
                U(nt).ons = data(1).ons;
                U(nt).dur = data(1).dur;
                U(nt).dt = data(1).dt; 
                microtime  = false;
            end
        end
    end
    
    U(1).grid = grid;
    U(1).gridpoints = argmatrix;
    U(1).grididx= idxmatrix;
    U(1).microtime = microtime;
    if ~microtime
        disp('microtime bins are not correctly set');
    end
    
    save(precomp_file,"U","grid","precomp_file");
    
    
    
end

end

