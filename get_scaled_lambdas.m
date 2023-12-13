function lambdas = get_scaled_lambdas(data_to_resample, mean_var)
% For an n_vars x n_timepoints dataset, compute varaince for each variable
% then rescale so that the variances across regions is distributed as
% Poi(mean_var).

norm_variances = zscore(var(data_to_resample, 0, 2));
lambdas = norm_variances * sqrt(mean_var) + mean_var;

end