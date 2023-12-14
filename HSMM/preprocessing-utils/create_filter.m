function f = create_filter(t_0, t, alpha, nv)
% Creates filter for filtered marked point process (fmpp).
% 
% PARAMETERS
% ----------
% t_0   : Reference timepoint/displacement.
% t     : n-vector of timepoints (note that these values are with respect
%         to 0, not t_0. 
% alpha : decay constant.
% Name-Value Pair (nv)
%   'renorm' : (1 | 0 (default)), specify whether or not to renormalize the
%              output so that it sums to 1.
% 
% RETURNS
% -------
% f : filter of length n.

arguments
    t_0
    t 
    alpha 
    nv.renorm = false
end

assert(isscalar(alpha), 'alpha must be a scalar.');
f = alpha * exp(alpha * (t - t_0));
f(t > t_0) = 0;

if nv.renorm, f = f ./ sum(f); end

end
