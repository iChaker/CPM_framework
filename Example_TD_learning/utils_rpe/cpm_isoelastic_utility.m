function out = cpm_isoelastic_utility(c, eta)
%ISOELASTIC_UTILITY, isoelastic utility function, for wealth (c) > 0 and real parameter
% eta. 
    if eta == 1
        out = log(c);
    else
        out = ((c.^(1 - eta) - 1 ) ./ (1 - eta));
    end
end