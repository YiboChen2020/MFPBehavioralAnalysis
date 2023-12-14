function B = get_overlap_mat(x, y, type)
% Get overlap between two state/bhv trajectories (pairwise between each
% state/behavior in each timeseries). 
%
% PARAMETERS
% ----------
% x    : k-vector, each of whose elements is a member of {r_1, r_2, ...,
%        r_n}, where r_i is a real number associated with the i_th behavior
%        or state (out of n such states/behaviors). 
% y    : k-vector, each of whose elements is a member of {q_1, q_2, ...,
%        q_m}, where q_j is a real number associated with the j_th behavior
%        or state (out of m such states/behaviors).
% type : {'cosine' | 'count'}. In either case, the function will compute
%        the overlap between the i_th and j_th state/behavior of x and y
%        respectively, for all i in {1, ..., n} and for all j in {1, ...,
%        m}. For each such pairing, the vectors x and y are converted into
%        Boolean valued vectors x* and y*, where the the p_th element of x*
%        is 1 if it equals r_i and 0 otherwise (and similarly, the p_th
%        element of y* is 1 if it equals q_j and 0 otherwise). If 'cosine',
%        the cosine similarity of x* and y* is then computed and stored as
%        the i_th,j_th element of the output matrix B. If 'count', let h be
%        the elementwise product of x* and y*. The i_th, j_th element of B
%        is then sum(h).
%
% RETURNS
% -------
% B : n x m matrix whose i_th,j_th element is the overlap between the i_th
%     bhv/state of x and the j_th bhv/state of y.
%
% Author: Jonathan Chien


if strcmp(type, 'cosine')
    x_list = unique(x);
    n_x = length(x_list);
    y_list = unique(y);
    n_y = length(y_list);
    
    B = nan(n_x, n_y);
    for i_x = 1:n_x
    for i_y = 1:n_y
        B(i_x, i_y) = cos_sim(double(x == x_list(i_x)), double(y == y_list(i_y)));
    end
    end
elseif strcmp(type, 'count')
    B = get_occupancy_mat(x, y);
else
    error("Unrecognized value for type.")
end

end
