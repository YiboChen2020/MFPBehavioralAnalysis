function expanded_labels = expand_bin_labels(labels, bin_width, step_size)
% For an n-vector whose i_th element is the latent state label of the i_th
% bin, return an (n*m)-vector, where m is the bin width (for step size also
% equal to m), where the i_th element is the latent state label of the i_th
% timepoint (assigned by propagating the bin label over all the individual
% timepoints within the bin).
%
% Author: Jonathan Chien. 11/1/22.

assert(bin_width == step_size, ...
       "Bin width must equal step size (i.e., no overlap of bins), else" + ...
       " it may be possible for a single timepoint to appear in multiple" + ...
       " bins, each potentially with a different label.")

expanded_labels = repelem(labels, bin_width);

end