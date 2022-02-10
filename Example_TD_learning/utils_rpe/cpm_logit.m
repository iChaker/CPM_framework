function out = cpm_logit(a)
%LOGIT Calculates the logit of a

% Added eps, to prevent over / underflow for values of 0 and 1.
out = log( (a + eps) ./ (1 - a + eps));

end
