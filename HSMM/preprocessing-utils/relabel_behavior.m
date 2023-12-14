function session_raw = relabel_behavior(session_raw, old, new)
% session_raw is a struct, not a cell or cell array containing struct(s).
% Author: Jonathan Chien

session_raw.behaviors = strrep(session_raw.behaviors, old, new);

end