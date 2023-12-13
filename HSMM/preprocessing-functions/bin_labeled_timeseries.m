function bin_labeled_timeseries(timeseries, labels, dim, bin_width, step_size)

arguments
    timeseries
    labels
    dim
    bin_width
    step_size
    nv.
end



% Assert nonoverlapping bins to ensure that transition between behaviors
% falls only in one bin.
assert(bin_width == step_size);

% Bin timeseries. 
binned_timeseries = bin_nd_data(timeseries, dim, bin_width, step_size, ...
                                'trim_from', 'back', 'bin_excess', 'excess_only');

% Determine bin index of all transition points.
trans_ind = find(diff(labels)) + 1;
trans_bin_ind = trans_ind ./ bin_width;
if any(mod(trans_bin_ind, 1) == 0)
    bin_starts_w_bhv = find(mod(trans_bin_ind, 1) == 0);
    disp('At least one bin overlaps with the start of a behavior.');
end
trans_bin_ind = ceil(trans_bin_ind);


binned_labels = nan(size(binned_timeseries, dim), 1);
for i_epoch = 1:length(trans_bin_ind)
    

end





