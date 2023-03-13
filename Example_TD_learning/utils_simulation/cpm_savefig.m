function cpm_savefig(fig, fname)

[~, ~, extension] = fileparts(fname);

if strcmp(extension, '.png')
    
    print(fig, fname,'-dpng', '-r700', '-painters');

elseif strcmp(extension, '.pdf')
    
    set(fig,'Units','points');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','points','PaperSize', [pos(3), pos(4)])
    print(fig, fname, '-dpdf', '-r600',  '-bestfit');
else
    error('%s - not yet implemented', extension)

end
   


end



