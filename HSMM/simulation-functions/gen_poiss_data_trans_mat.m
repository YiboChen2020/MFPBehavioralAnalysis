function [sim_data, states, degenerate, seed] ...
    = gen_poiss_data_trans_mat(activity_mat, trans_mat, n_timepoints, nv)
%
% Uses an empirical transition matrix and mean activity matrix to generate
% a timeseries of counts based on draws from a Poisson distribution.
%
% PARAMETERS
% ----------
% activity_mat : k x m matrix whose i_th, j_th element is the mean activity
%                of the m_th region in the k_th state/behavior etc.
% trans_mat    : k x k transition matrix whose i_th, j_th element is the
%                probability of transitioning from the i_th state to the
%                j_th state.
% n_timepoints :
% Name-Value Pair (nv)
%   'starting_bhv_idx'  : ([] (default) | scalar integer) Specify index of
%                         behavior to start from. Else, select a behavior 
%                         at random.
%   'add_gaussian_nosie : (false (default) | 4-cell-array) Specify a 4
%                         element cell array whose
%                         cells contain, respectively, the variance of
%                         Gaussian noise to be overlayed over the simulated
%                         data, the Poisson parameter for resampling, the
%                         seed, and a Boolean value indicating whether to
%                         pool all regions and resample (if true), or
%                         resample each region separately.
%
% RETURNS
% -------
% sim_data : m x n array of simulated data.
% states   : n-vector whose i_th element is the state label of the i_th
%            column of states.
%
% Author: Jonathan Chien. 10/15/22.

arguments
    activity_mat
    trans_mat
    n_timepoints
    nv.starting_bhv_idx = []
    nv.add_gaussian_noise = false % If true, specify 2-vector where first element is gaussian variance, and second is Poisson parameter for resampling, e.g., [5 10]
    nv.remove_degenerate = false
    nv.seed = 'default'
end

% Check if any diagonal elements of transition matrix are 1 and warn if so.
degenerate = ismembertol(diag(trans_mat), 1, 'Datascale', 1);
if any(degenerate)
    if nv.remove_degenerate
        warning("State(s) %s has/have self-transition probability of 1 and " + ...
                "will be removed.", ...
                num2str(find(ismembertol(diag(trans_mat), 1))));

        % Remove degenerate state(s) from transition matrix and renormalize.
        trans_mat(degenerate, :) = [];
        trans_mat(:, degenerate) = [];
        trans_mat = trans_mat ./ sum(trans_mat, 2);

        % Remove degenerate state(s) from activity matrix.
        activity_mat(degenerate, :) = [];
    else
        warning("State(s) %s has/have self-transition probability of 1 but will not be removed. " + ...
                "The model will be unable to leave this/these state(s) if it " + ...
                "enters it/them.", num2str(find(ismembertol(diag(trans_mat), 1))));
    end
end

% Preallocate.
[n_bhvs, n_regions] = size(activity_mat);
sim_data = NaN(n_regions, n_timepoints);
states = NaN(n_timepoints, 1);

% Set starting behavior randomly if not specified by user.
if ~isempty(nv.starting_bhv_idx)
    i_bhv = nv.starting_bhv_idx;
else
    i_bhv = randi(n_bhvs, 1, 1);
end

% Set seed.
rng(nv.seed);
seed = rng;

% Generate data.
for i_timepoint = 1:n_timepoints
    % Draw data for current timestep.
    sim_data(:,i_timepoint) = poissrnd(activity_mat(i_bhv,:), 1, n_regions);
    states(i_timepoint) = i_bhv;

    % Draw bhv for next step based on transition matrix.
    i_bhv = find(mnrnd(1, trans_mat(i_bhv,:)));
end

% Optionally corrupt with Gaussian noise and then resample back to Poisson.
if iscell(nv.add_gaussian_noise)
    sim_data = corrupt_and_resample(sim_data, ...
                                    nv.add_gaussian_noise{1}, ...
                                    nv.add_gaussian_noise{2}, ...
                                    nv.add_gaussian_noise{3}, ...
                                    nv.add_gaussian_noise{4});
end

end
