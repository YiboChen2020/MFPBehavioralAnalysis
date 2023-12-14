function occupancy = get_occupancy_mat(x, y)
% Gets count based overlap between two state/bhv trajectories (pairwise
% between each state/behavior in each timeseries).
%
% PARAMETERS
% ----------
% x    : k-vector, each of whose elements is a member of {r_1, r_2, ...,
%        r_n}, where r_i is a real number associated with the i_th behavior
%        or state (out of n such states/behaviors). 
% y    : k-vector, each of whose elements is a member of {q_1, q_2, ...,
%        q_m}, where q_j is a real number associated with the j_th behavior
%        or state (out of m such states/behaviors).
% type : The function will compute
%        the overlap between the i_th and j_th state/behavior of x and y
%        respectively, for all i in {1, ..., n} and for all j in {1, ...,
%        m}. For each such pairing, the vectors x and y are converted into
%        Boolean valued vectors x* and y*, where the the p_th element of x*
%        is 1 if it equals r_i and 0 otherwise (and similarly, the p_th
%        element of y* is 1 if it equals q_j and 0 otherwise). Let h be
%        the elementwise product of x* and y*. The overlap between the i_th
%        and j_th states/bhvs is then sum(h) and is stored as the i_th,
%        j_th element of B.
%
% RETURNS
% -------
% B : n x m matrix whose i_th,j_th element is the overlap between the i_th
%     bhv/state of x and the j_th bhv/state of y.
%
% Author: Jonathan Chien

% These lists will be returned in sorted order by default.
x_list = unique(x);
y_list = unique(y);

% Count overlap between each pairing of states/bhvs in x and y.
occupancy = nan(length(x_list), length(y_list));
for i_x = 1:length(x_list)
for i_y = 1:length(y_list)
    occupancy(i_x, i_y) ...
        = sum(ensure_col(x) == x_list(i_x) & ensure_col(y) == y_list(i_y));
end
end

end
