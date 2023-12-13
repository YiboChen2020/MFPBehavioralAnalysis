function f = create_filter(t_0, t, alpha, nv)
% Shifted exponential distribution reflected across y-axis.
% 
% PARAMETERS
% ----------
% t_0   : Reference timepoint, equal to the value by which we displace the
%         reflected exponential distribution.
% t     : n-vector of timepoints (values in the support of the reflected
%         exponential distribution; note that these values are with respect
%         to 0, not t_0.
% alpha : decay constant.
% Name-Value Pair (nv)
%   'renorm' : (1 | 0 (default)), specify whether or not to renormalize the
%              n draws so that they sum to 1 (treating them more like
%              probability masses).
% 
% RETURNS
% -------
% f : n-vector consisting of n draws from reflected exponential
%     distribution (possibly renormalized so that discrete values sum to
%     1).

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

% From the 2019 paper, but this form seems to be based on the CDF so not
% exactly sure what's going on here:
% f = (1 - (exp(-alpha * (t - t_0)))) / (exp(-alpha * t_0) / alpha);
% f(t > t_0) = 0;
