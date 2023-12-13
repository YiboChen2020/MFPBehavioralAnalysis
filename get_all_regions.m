function all_regions = get_all_regions(path_to_data, ind_range)
% Search all sessions from all subjects and get a list of all regions.
%
% PARAMETERS
% ----------
% path_to_data : Path to directory containing data.
% ind_range    : 2-vector specifying range of subject indices to try.
%
% RETURNS
% -------
% all_behaviors : n x 1 cell array, where n = the number of unique regions
%                 identified, and the i_th cell contains the name of the
%                 i_th region.

addpath(path_to_data)

% Initialize list variable.
all_regions = [];

% Search across all subjects and sessions.
for i_subject = ind_range(1):ind_range(2)
    % Load subject's data.
    try
        data = load(sprintf('MFP-ERa-GC%d.mat', i_subject));
    catch
        continue
    end

    all_regions = horzcat(all_regions, data.PETH.regions);
end

all_regions = unique(all_regions);

end