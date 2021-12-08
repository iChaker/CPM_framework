function U = cpm_precompute(model,grid,fixedparams,data,precomp_file)
%CPM_PRECOMPUTE Summary of this function goes here
%   Detailed explanation goes here



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

    U(length(ons)) = struct();
    S= {};
    S.subs = repmat({':'},1,ndims(precomputed));
    S.type = '()';

    for nt = 1 : length(ons)
        S.subs{1}=nt;
        signals=squeeze(subsref(precomputed,S));
        
        signals1D = zeros(length(idxmatrix),1);
        for i=1:length(idxmatrix)
            idx = num2cell(idxmatrix(i,:));
            signals1D(i)=signals(idx{:});
        end
        
    
        U(nt).signals = signals;
        U(nt).signals1D = signals1D;
        U(nt).ons = data.ons(nt);
        U(nt).dur = data.dur(nt);
        U(nt).dt = data.dt(nt); 
    end
    
    U(1).grid = grid;
    U(1).gridpoints = argmatrix;
    U(1).grididx= idxmatrix;
    save(precomp_file,"U","grid","precomp_file");
    
    
    
end

end

