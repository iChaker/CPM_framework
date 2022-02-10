function out = cpm_logit_inv(a)
%LOGIT_INV Summary calculates the inverse logit of a
out = exp(a) ./ (1 + exp(a));
end
