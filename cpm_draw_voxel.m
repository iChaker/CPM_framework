function z = cpm_draw_voxel(PRF,voxel,dim1,dim2,fixed_values)
%CPM_DRAW_VOXEL Summary of this function goes here
%   Detailed explanation goes here
%   Voxel:    voxel index
%   dim1 - dim2 :     parameters names to show i.e 'tau' 'eta'
%   fixed_values: array of fixed values for the rest of the parameters
%   (defaults to array of zeros)
%   Make sure dim1 and dim2 are ordered ('tau' 'eta' is correct 'eta' 'tau'
%   gives a wrong graph)

grid = PRF.U(1).grid;

d1 = grid.(dim1);
d2 = grid.(dim2);


b=80;


% Make correct coords

p_names = fieldnames(grid);
j=1;

try 
    fixed_values;
catch
    fixed_values= zeros(1,length(p_names)-2);
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

figure;

z  = feval(PRF.M.IS, PRF.Ep{voxel}, PRF.M, PRF.U, 'get_response', coords);
z  = reshape(z,b,b);

ticks = linspace(1,b,17);

colormap(jet);
imagesc(z);
%[ 'voxel:' num2str(voxel) 'at' num2str(fixed_values)]
title([ 'Voxel: ' num2str(voxel) ' slice for: ' num2str(fixed_values) ] ,'FontSize',16);
set(gca,'YDir','normal','XTick',ticks,'YTick',ticks);
set(gca,'XTickLabel',linspace(d1(1),d1(2),length(ticks)),'YTickLabel',linspace(d2(1),d2(2),length(ticks))); 
xlabel(dim1);
ylabel(dim2);


end

