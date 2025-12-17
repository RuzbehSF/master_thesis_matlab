
%% IRI_OF_LTPP
% -------------------------------------------------------------------------
% PURPOSE:
%   Import IRI results already computed in ProVAL for the LTPP profile used
%   in Section 1. No IRI computation here—just parsing, unit sanity, save.
%
% OUTPUTS:
%   ./00_Outputs/05_IRI_Analysis/TrueProfile/iri_true_profile.mat
%   ./00_Outputs/05_IRI_Analysis/TrueProfile/iri_true_profile.csv
% -------------------------------------------------------------------------

close all; clc;

%% 0) Config — set your ProVAL IRI CSV here
iri_cfg = struct( ...
    'csv_path', "LTPP_Input/IRI/profile_9_IRI.csv", ...   % ProVAL Ride Quality export
    'lane',     "Center" ...                               % lane/series label to prefer
);

%% 1) Read CSV (preserve headers exactly as written)
tbl = readtable(iri_cfg.csv_path, 'VariableNamingRule','preserve');

% --- Column detection (robust to ProVAL header wording) ---
% Start distance
startVar = pickVar(tbl, [ ...
    "Start Distance (m)","Start (m)","Begin (m)","From (m)","Start","Begin","From" ...
]);

% End/Stop distance (ProVAL often uses "Stop Distance (m)")
endVar = pickVar(tbl, [ ...
    "Stop Distance (m)","End Distance (m)","End (m)","Stop (m)","To (m)","End","Stop","To","Finish" ...
]);

% Length (as a fallback to build End = Start + Length)
lenVar = pickVar(tbl, [ ...
    "Length (m)","Segment Length (m)","Length","Segment Length" ...
]);

% IRI column — prefer lane-specific first, then generic
lane = string(iri_cfg.lane);
iriVar = pickVar(tbl, [ ...
    lane + " IRI (m/km)", lane + " IRI (in/mi)", lane + " IRI", ...
    lane + " Elevation - IRI (m/km)", lane + " Elevation - IRI", ...
    "Center Elevation - IRI (m/km)", "Center Elevation - IRI", ...
    "IRI " + lane + " (m/km)", "IRI " + lane + " (in/mi)", "IRI " + lane, ...
    "IRI (m/km)","IRI (in/mi)","IRI" ...
]);

if strlength(startVar) == 0
    error('iri_of_LTPP:StartColumnNotFound', ...
        'Could not find a Start distance column in %s', iri_cfg.csv_path);
end
if strlength(iriVar) == 0
    error('iri_of_LTPP:IRIColumnNotFound', ...
        'Could not find an IRI column in %s', iri_cfg.csv_path);
end
if (strlength(endVar) == 0) && (strlength(lenVar) == 0)
    error('iri_of_LTPP:EndOrLengthNotFound', ...
        'Need either an End/Stop distance column or a Length column; none found in %s', iri_cfg.csv_path);
end

%% 2) Pull and clean vectors
x0_raw = tbl.(startVar);

if strlength(endVar) > 0
    x1_raw = tbl.(endVar);
else
    L_raw  = tbl.(lenVar);
    x1_raw = x0_raw + L_raw;   % infer end from start + length
end

iri_raw = tbl.(iriVar);

good  = isfinite(x0_raw) & isfinite(x1_raw) & isfinite(iri_raw);
x0_raw = x0_raw(good);
x1_raw = x1_raw(good);
iri_raw= iri_raw(good);

% --- Units ---
% Start/End to meters (read from header text "(m)", "(ft)", "(mi)")
u_start = headerUnits(tbl, startVar);
if strlength(endVar) > 0
    u_end = headerUnits(tbl, endVar);
else
    % If no endVar, we formed x1 from start + length; use length column units
    u_end = headerUnits(tbl, lenVar);
end

x0 = convert_distance_to_m(x0_raw, u_start);
x1 = convert_distance_to_m(x1_raw, u_end);

% IRI to m/km (handle in/mi)
u_iri = headerUnits(tbl, iriVar);
if contains(lower(u_iri), "in/mi")
    iri_m_per_km = iri_raw * (0.0254 / 1609.344) * 1000;  % in/mi → m/km
else
    iri_m_per_km = iri_raw;  % assume already m/km
end

% Keep positive-length segments
len_m = x1 - x0;
pos   = isfinite(len_m) & (len_m > 0);
x0    = x0(pos);
x1    = x1(pos);
len_m = len_m(pos);
iri_m_per_km = iri_m_per_km(pos);

if isempty(x0)
    error('iri_of_LTPP:NoSegments', 'No valid IRI segments after cleaning (nonpositive or NaN lengths).');
end

%% 3) Length-weighted overall IRI
overall_m_per_km = sum(iri_m_per_km .* len_m) / sum(len_m);

%% 4) Save artifacts where your comparer expects them
outDir = fullfile('00_Outputs','05_IRI_Analysis','TrueProfile');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

T = table(x0, x1, len_m, iri_m_per_km, ...
    'VariableNames', {'Start_m','End_m','Length_m','IRI_m_per_km'});

writetable(T, fullfile(outDir,'iri_true_profile.csv'));

IRI = struct();
IRI.meta = struct('source','ProVAL/LTPP','csv',string(iri_cfg.csv_path),'lane',lane);
IRI.overall_m_per_km = overall_m_per_km;
IRI.segment_length_m = median(len_m,'omitnan');
IRI.segments = T;

save(fullfile(outDir,'iri_true_profile.mat'), 'IRI', '-v7.3');

fprintf('[IRI LTPP] Overall IRI = %.3f m/km | segments: %d | saved → %s\n', ...
    overall_m_per_km, height(T), outDir);

%% ---------------- local helpers ----------------
function name = pickVar(tbl, candidates)
% Return the first matching variable name (exact, then contains; case-insensitive)
    if isstring(candidates) || ischar(candidates)
        candidates = string(candidates);
    end
    vars = string(tbl.Properties.VariableNames);

    % 1) exact match
    for c = candidates
        hit = vars(strcmpi(vars, c));
        if ~isempty(hit)
            name = hit(1);
            return;
        end
    end
    % 2) contains
    low = lower(vars);
    for c = candidates
        if strlength(c) == 0
            continue;
        end
        w = lower(c);
        hit = vars(contains(low, w));
        if ~isempty(hit)
            name = hit(1);
            return;
        end
    end
    name = "";
end

function u = headerUnits(tbl, varName)
% Extract "(units)" either from VariableDescriptions or from the header text
    txt = "";
    vdesc = tbl.Properties.VariableDescriptions;
    if ~isempty(vdesc)
        idx = find(strcmp(tbl.Properties.VariableNames, varName), 1);
        if ~isempty(idx) && ~isempty(vdesc{idx})
            txt = string(vdesc{idx});
        end
    end
    if txt == ""
        txt = string(varName);
    end
    tok = regexp(txt, "\(([^)]+)\)", "tokens", "once");
    if isempty(tok)
        u = "";
    else
        u = string(tok{1});
    end
end

function x = convert_distance_to_m(x, u)
% Convert distance vector x to meters based on a unit string u
    switch lower(string(u))
        case "m"   % already meters
        case "km", x = x*1000;
        case "ft", x = x*0.3048;
        case "mi", x = x*1609.344;
        otherwise  % unknown → assume meters
    end
end