function remove_ind = bhv_name_to_remove_ind(bhvs_to_exclude, ...
                                             bhv_name_codes, ...
                                             unique_bhvs)
% ------------------------------------- 
% Warning: this function is deprecated. 
% ------------------------------------- 
%
% Routine/subroutine to convert string names of behaviors to be excluded to
% column indices of counts matrix/row indices of activity matrix (note that
% column/row indices need not match up with the number code of the
% behavior).
%
% Note for pooling sessions for one subject: it is possible to include
% behaviors that are not actually present in a given session's data (they
% will simply be ignored). For example, if we call labeled_timeseries for
% all sessions for a given subject, all behaviors will be pooled across all
% sessions first, and a unique number then assigned to each unique
% behavior. Thus it may be the case that a particular session (e.g.,
% session 1) does not feature a given behavior (e.g., Intro_Toy, with
% corresponding code 5), since Intro_Toy occurred in a different session.
% The name "Intro_Toy" and the corresponding code, 5, can still be passed
% into this function, however, and it will simply be ignored. This is
% useful, as we can first get the list of all behaviors that occured across
% all sessions for a subject, define some subset of the behaviors (e.g.,
% Investigate, Attempted Attack, and Attack), and then get the set
% difference of the full list and this subset. The set difference is the
% list of all behaviors to be excluded. Not all of these behaviors to be
% excluded will have occurred in all of the sessions, but we can still pass
% in the same list of behaviors to be excluded for all sessions.
%
% PARAMETERS
% ----------
% bhvs_to_exclude : Cell array whose i_th element contains the name of one
%                   of the behaviors that we would like to exclude. Note
%                   that it is possible to include behaviors that are not
%                   actually present in a given session's data (they will
%                   simply be ignored). 
% bhv_name_codes  : Vector whose i_th element is the number code associated
%                   with the i_th element of bhvs_to_exclude. As with
%                   bhvs_to_exclude, extra codes will be ignored.
% unique_bhvs     : Vector of number codes that actually occured in the
%                   data being processed.
%
% RETURNS
% -------
% remove_ind : Indices of behaviors to be removed.
%
% Author: Jonathan Chien

if ~iscell(bhvs_to_exclude), remove_ind = []; return; end

% Get code numbers of behaviors to be removed. Note that the
% bhv_name_codes field has as many names (in order by code number
% and by alphetical order) as codes. Thus we can work directly with
% the names here.
codes_to_exclude ...
    = find(ismember(bhv_name_codes(:,1), bhvs_to_exclude));

% Get index of count matrix from number code. Note that unique_bhvs are
% sorted when they are derived in the parent function.
if isrow(codes_to_exclude), codes_to_exclude = codes_to_exclude'; end
if iscolumn(unique_bhvs), unique_bhvs = unique_bhvs'; end
remove_ind = find(sum(sort(codes_to_exclude) == unique_bhvs, 1));

end
