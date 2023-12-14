function [sessions_raw, selected_session_ids] ...
    = extract_single_subj(path_name, i_subject, varargin)
% Extracts all raw data from all sessions (or subset of sessions) for one
% subject. Optional visualization of raw traces.
%
% PARAMETERS
% ----------
% path      : String, path to the directory containing raw data.
% i_subject : Scalar integer, index of the subject to be analyzed.
% Name-Value Pairs (nv)
%   'session'     : ('all' (default) | string | vector), set as 'all' to
%                   load data from all sessions. Or, pass in the desired
%                   session ID as a string, or pass in a vector of indices
%                   of sessions we wish to load; indices must be from the
%                   set {1:n}, where n = the number of sessions for the
%                   current subject. The length of this vector is equal to
%                   n_sessions_selected (see RETURNS).
%   'region_names' : (1 | 0 (default)), if region names are on the same path
%                   as the neural data, set true to load them and add them
%                   to the returned sessions_raw struct.
%   'remove_nan'  : (1 (default) | 0), specify whether to remove regions
%                   that have NaN values.
%   'nan_warning' : (1 (default) | 0), specify whether to issue warning
%                   if regions were removed due to NaN values. Ignored if
%                   'remove_nan' = false.
%   'zscore'      : (1 | 0 (default)), specify whether or not to z-score
%                   each region's timeseries separately.
%   'smooth'      : (5 | scalar | 0 (default)), specify the window size for
%                   gaussian smoothing, or set to 0/false to suppress
%                   smoothing.
%   'plot'        : (1 | 0 (default)), specify whether or not to visualize
%                   raw traces.
%   'figure'      : (1 (default) | 0), specify whether or not to plot into
%                   new figure. Ignored if 'plot' is false.
%   'plot_scale'  : Scalar value multiplying the raw traces; can be used
%                   to scale down traces to prevent overlap. Ignored if
%                   'plot' is false.
%   'window'      : 2-vector whose 1st and 2nd elements are the indices
%                   of the 1st and 2nd timepoints, respectively, of the
%                   window we would like to visualize. Ignored if 'plot' is
%                   false.
%
% RETURNS
% -------
% sessions_raw         : n_sessions_selected x 1 cell array. The k_th cell
%                        contains a scalar struct with the following
%                        fields:
%   .Lfold                : Recorded GCaMP traces
%   .t                    : Time stamp of each GCaMP imaging frame
%   .behaviors            : Annotated behavior
%   .session_id           : E.g., '200801'
%   .Fstart               : Starting frame of behaviors
%   .FStop                : Stopping frame of behaviors
%   .regions              : Brain region
%   .FL                   : Sampling rate, 25 Hz --> 40 ms
%   .n_regions_retained   : n_sessions_selected x 1 array whose i_th
%                           element is the number of regions retained after
%                           possible NaN removal for the i_th session.
%   .nan_regions          : This is an n_sessions_selected x 1 cell array
%                           whose i_th cell contains a vector whose
%                           elements are the indices of regions removed due
%                           to containing NaN values. If 'nan_regions' =
%                           false, all cells will be empty.
% selected_session_ids : n_sessions_selected x 1 cell array whose i_th cell
%                        contains the session ID of the i_th session
%                        selected/to be returned.
%
% Author: Jonathan Chien. 5/19/22. 


%% Parse inputs

p = inputParser;
addRequired(p, 'path', @ischar);
addRequired(p, 'i_subject', @(x) isscalar(x) && all(floor(x) == x));
addParameter(p, 'session', 'all', @(x) ischar(x) || isstring(x) || iscell(x) || all(floor(x) == x))
addParameter(p, 'combine_sessions', false)
addParameter(p, 'bhvs_to_exclude', [])
addParameter(p, 'region_names', false, @islogical)
addParameter(p, 'remove_nan', true, @islogical)
addParameter(p, 'nan_warning', true, @islogical)
addParameter(p, 'zscore', false)
addParameter(p, 'smooth', false)
addParameter(p, 'plot', false, @(x) islogical(x))
addParameter(p, 'figure', true, @(x) islogical(x))
addParameter(p, 'plot_scale', 0.3, @(x) isscalar(x))
addParameter(p, 'window', [1000 3000], @(x) length(x) == 2 && all(floor(x) == x))
parse(p, path_name, i_subject, varargin{:});
nv = p.Results;


%% Load data

% Add directory containing data to path.
c = pathsep;
path_str = [c path c];
on_path = contains(path_str, [c path_name c], 'IgnoreCase', ispc);
if ~on_path, addpath(path_name); end

% Set subject index and load data. Since there are technically many ways to
% generate an error with load, this setup here probably isn't a great way
% to do things, but we will treat any errors as due to the filename not
% existing and throw an exception in invoking function named (as of 6/1/22)
% snr_all_subjects_session; this will at least cause errors within this
% function not related to file-loading to be caught.
filename = sprintf('MFP-ERa-GC%d.mat', i_subject);
try
    data_as_loaded = load(filename);
catch
    exception = MException("extract_single_subject:file_not_found", ...
                           "Filename %s not found", filename);
    throwAsCaller(exception);
end

% Optionally load region names if available.
if nv.region_names
    try
        region_names = load("region-names.mat");
        region_names = region_names.regions;
    catch
        warning("'region_name' is true, but there might not be region " + ...
                "names available.")
    end
end

% Determine number of sessions and get session IDs.
fnames = fieldnames(data_as_loaded.Raw);
n_sessions = length(fnames) - 6; % There are six fields before the sessions
all_session_ids = fnames(7 : end);


%% Select data from all sessions or from subset of sessions

% Set the indices of the sessions we would like to run. Set 'all' or empty
% array to load all sessions. Or can specify session IDs, or session indices.
if iscell(nv.session) % Session ID
    % Ensure cell array contains strings/characters of session IDs.
    assert(all(cellfun(@isstring, nv.session)) || all(cellfun(@ischar, nv.session)), ...
           "If passing in cell array of session IDs, must be all strings or characters."); 
    
    % Get index of each session.
    session_ind = NaN(length(nv.session), 1);
    for i_session = 1:length(session_ind)
        session_ind(i_session) ...
            = find(~cellfun(@isempty, (regexp(all_session_ids, nv.session{i_session}))));
    end
    % Mutate nv.session into numeric array of session indices.
    nv.session = session_ind;
elseif strcmp(nv.session, 'all') || isempty(nv.session)
    nv.session = 1:n_sessions; 
else % Session indices
    assert(~isempty(nv.session))
    assert(isnumeric(nv.session))
    assert(all(floor(nv.session) == nv.session)) % Elements must be integer indices
end
n_selected_sessions = length(nv.session);

% Preallocate across selected sessions.
sessions_raw = cell(n_selected_sessions, 1);

% Load data from selected sessions. 
for i_session = 1:n_selected_sessions
    sessions_raw{i_session} ...
        = data_as_loaded.Raw.(all_session_ids{nv.session(i_session)});
    sessions_raw{i_session}.session_id = {all_session_ids{nv.session(i_session)}(2:end)};
    if nv.region_names
        sessions_raw{i_session}.region_names = region_names;
    end
end

% Get IDs of selected sessions.
selected_session_ids = all_session_ids(nv.session);


%% Remove any behaviors specified for removal

for i_sess = 1:n_selected_sessions
    % Get indices of markers for epochs to be removed.
    remove_ind = [];
    for i_bhv = 1:length(nv.bhvs_to_exclude)
        remove_ind = [remove_ind; find(strcmp(nv.bhvs_to_exclude, sessions_raw{i_sess}.behaviors))];
    end

    % Remove markers. 
    sessions_raw{i_sess}.behaviors(remove_ind) = [];
    sessions_raw{i_sess}.Fstart(remove_ind,:) = [];
    sessions_raw{i_sess}.Fstop(remove_ind,:) = [];
end


%% Check for NaNs in Lfold data

% Preallocate array whose i_th cell contains indices of regions that were
% removed for the i_th session.
nan_regions = cell(n_selected_sessions, 1);

% Preallocate array whose i_th cell contains the number of regions retained
% for the i_th session, after potentially dropping any regions.
n_regions_retained = NaN(n_selected_sessions, 1);

% Check all regions for each session, and remove regions with NaN values if
% desired. Also remove the corresponding region name. Optionally issue
% warning if regions removed.
for i_session = 1:n_selected_sessions
    nan_regions{i_session} = find(any(isnan(sessions_raw{i_session}.Lfold), 2));
    if nv.remove_nan && any(nan_regions{i_session})
        sessions_raw{i_session}.Lfold(nan_regions{i_session},:) = [];
        if nv.region_names, sessions_raw{i_session}.region_names(nan_regions{i_session}) = []; end
        if nv.nan_warning
            warning('Region(s) %d from session %s have NaN values and will be removed.', ...
                    nan_regions{i_session}, selected_session_ids{i_session})
        end
    end
    
    % Note: in the future, it is possible that there are other arrays that
    % will need to be similarly modified.
    
    % Record number of regions retained for this session.
    n_regions_retained(i_session) = size(sessions_raw{i_session}.Lfold, 1);

    % Add as field of sessions_raw.
    sessions_raw{i_session}.nan_regions = nan_regions{i_session};
    sessions_raw{i_session}.n_regions_retained = n_regions_retained(i_session);
end


%% Optionally z-score and smooth data

if nv.zscore
    for i_session = 1:n_selected_sessions
        sessions_raw{i_session}.Lfold = zscore(sessions_raw{i_session}.Lfold, 0, 2);
    end
end

if nv.smooth
    for i_session = 1:n_selected_sessions
        sessions_raw{i_session}.Lfold ...
            = smoothdata(sessions_raw{i_session}.Lfold, 2, 'gaussian', nv.smooth);
    end
end


%% Optionally combine all sessions into one
% Note that indivual values (e.g. sampling rate) will also be concatenated,
% such that if e.g. three sessions are combined, there will be a three
% element vector containing the sampling rate from each session.

if nv.combine_sessions && n_selected_sessions > 1
    % Copy the first session, then add the others.
    combined_sessions_raw = sessions_raw{1};

    % Get fieldnames (everything that needs to be concatenated).
    sessions_raw_fnames = fieldnames(sessions_raw{1});

    % 5/23/23: After Dayu's group modified/corrected some of the errors in
    % the data related to sampling, there appeared new fields
    % Lfold_original and FL_original. These are to be excluded.
    sessions_raw_fnames = setdiff(sessions_raw_fnames, {'Lfold_original', 'FL_original'});
    
    for i_field = 1:length(sessions_raw_fnames)
        % Determine dimension along which to concatenate.
        if ismember(sessions_raw_fnames{i_field}, {'Lfold', 't'})
            cat_dim = 2;
        else
            cat_dim = 1;
        end

        for i_session = 2:n_selected_sessions
            % XXX: this will currently generate an error during attempt to
            % concatenate region names cell arrays along first dim if the
            % number of regions retained was different between sessions,
            % though the number of regions generally appears to be the same
            % across sessions.
            try
                combined_sessions_raw.(sessions_raw_fnames{i_field}) ...
                    = cat(cat_dim, ...
                          combined_sessions_raw.(sessions_raw_fnames{i_field}), ...
                          sessions_raw{i_session}.(sessions_raw_fnames{i_field}));
            catch
                keyboard
                % TODO: pad with nans or other placeholder to allow
                % concatenation.
            end
        end
    end

    % Reassigned combined struct as new sessions_raw. NB:
    % n_selected_sessions is also mutated to 1 here.
    sessions_raw = {combined_sessions_raw};
    n_selected_sessions = 1;
end


%% Optionally plot Lfold from all selected sessions for current subject

if nv.plot      
    if nv.figure, figure; end

    % Plot all selected sessions for current subject.
    for i_session = 1:n_selected_sessions
    for i_region = 1:n_regions_retained(i_session)   
       subplot(n_selected_sessions, 1, i_session)
       hold on
       
       plot(i_region + sessions_raw{i_session}.Lfold(i_region,:)/nv.plot_scale, ...
            'linewidth', 2);   
       
       xlim(nv.window)
       xlabel('Time')
       ylabel('Region index')
       title(sprintf('Session %s', selected_session_ids{i_session}))     
    end
    end; clear i_session i_region

    % Resize figure window.
    set(gcf, 'Position', [100, 500, 1250, n_selected_sessions*400]);   
end

end
