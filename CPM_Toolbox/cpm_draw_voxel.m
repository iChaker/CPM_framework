function z = cpm_draw_voxel(PRF,Pvoxel,dim1,dim2,fixed_values,bins,tile)
%cpm_draw_voxel draws population field over a slice of the parameter space
%   INPUT:
%   Voxel:    voxel index
%   dim1 - dim2 :     parameters names to show i.e 'tau' 'eta'
%   fixed_values: array of fixed values for the rest of the parameters
%   (defaults to array of zeros)
%   Make sure dim1 and dim2 are ordered ('tau' 'eta' is correct 'eta' 'tau'
%   gives a wrong graph)
%   bins: defines granularity of the image
%
%   OUTPUT:
%   z: weights of population field for the given slice and precision.

grid = PRF.U(1).grid;

d1 = grid.(dim1);
d2 = grid.(dim2);


b=bins;

% Make correct coords

p_names = fieldnames(grid);
j=1;

try 
    if isempty(fixed_values)
        fixed_values= zeros(1,length(p_names)-2);
    end
catch
    fixed_values= zeros(1,length(p_names)-2);
end

try 
    P=PRF.Ep{Pvoxel};
catch
    P= Pvoxel;
end



x_bins = linspace(d1(1),d1(2),b);
y_bins = linspace(d2(1),d2(2),b);

[x2,y2] = meshgrid(x_bins,y_bins);
xy=[x2(:) y2(:)];
xy=xy';

coords = zeros(length(xy),length(p_names));

for i=1:length(p_names)

        if isequal(p_names{i},dim1)
           coords(:,i) = xy(1,:);
        
        elseif isequal(p_names{i},dim2)
            coords(:,i) = xy(2,:);
        
        else
            coords(:,i)=fixed_values(j) * ones(length(xy),1);
            j=j+1;
        end

end

try tile;
   nexttile(tile); 
   c=9;
catch
   figure;
   c=17;
end
z  = feval(PRF.M.IS, 	P , PRF.M, PRF.U, 'get_response', coords);
z  = reshape(z,b,b);

ticks = linspace(1,b,c);

contourf(peaks)
colormap(jet);
imagesc(z);
%[ 'voxel:' num2str(voxel) 'at' num2str(fixed_values)]
if c==17
    try 
        title([ 'Voxel: ' num2str(Pvoxel) ' slice for: ' num2str(fixed_values) ] ,'FontSize',16); 
    catch
        title([ 'Gnerative Voxel, slice for: ' num2str(fixed_values) ] ,'FontSize',16);     
    end
    xlabel(dim1);
    ylabel(dim2);
end


set(gca,'YDir','normal','XTick',ticks,'YTick',ticks);
set(gca,'XTickLabel',linspace(d1(1),d1(2),length(ticks)),'YTickLabel',linspace(d2(1),d2(2),length(ticks))); 



end

