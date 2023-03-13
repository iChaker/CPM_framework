function output = cpm_grid_SAMSRF(freeparams,fixedparams,data)
% The model is now simply a look up...

appertures = data.appertures;
resolution = fixedparams.resolution; % because 0, 0 is the middle, we need to add to 
offset = ceil(resolution / 2);

signal = squeeze(appertures(:, freeparams.x + offset, freeparams.y + offset));

abs_max = max(abs(signal(:))) + eps; % if everything is 0
signal = signal(:) ./ abs_max;  

output=signal;

end




