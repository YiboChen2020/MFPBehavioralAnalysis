function sim_data = corrupt_and_resample(sim_data, noise_var, lambda, seed, pool_flag)
% Accepts an n_regions x n_timepoints array of simulated data, adds
% Gaussian noise, and then resamples back to Poisson.
% 
% PARAMETERS
% ----------
% sim_data  : m x n array of simulated data, m = n_regions, n =
%             n_timepoints.
% noise_var : Scalar factor that multiplies identity matrix to form
%             covariance matrix of isotropic Gaussian noise overlaid on
%             simulated data.
% lambda    : Poisson parameter for resamping after noise corruption.
% seed      : seed for resampling.
% pool_flag : {1 | 0}, specify whether to pool regions together and
%             resample across all regions (if true), or not (if false).
%
% Author: Jonathan Chien

% Add Gaussian noise.
sim_data = sim_data + randn(size(sim_data)) * sqrt(noise_var);

% Resample. Either pool data from all regions and resample together, or
% resample each region separately.
if pool_flag
    sim_data = resample2poisson(sim_data, lambda, seed);
elseif ~pool_flag
    for i_region = 1:size(sim_data, 1)
        sim_data(i_region,:) = resample2poisson(sim_data(i_region,:),lambda, seed);
    end
end

end