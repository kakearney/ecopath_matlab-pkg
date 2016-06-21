function A = unitconvert(A, oldunit, newunit, varargin)


% Parse and check input

p = inputParser;
p.addParameter('wwCfrac', 0.05,    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('c2n',     16/106,  @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('cmw',     12.0107, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('depth',   NaN,     @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.parse(varargin{:});

Opt = p.Results;

%------------------------------
% Calculate conversion factors
%------------------------------

% Table of ratios

mass = {...
    'tons wet weight'   1
    't ww'              1
    't C'               Opt.wwCfrac
    'g C'               Opt.wwCfrac * 1e6
    'mmol N'            (Opt.wwCfrac * Opt.c2n /Opt.cmw) * 1e9
    'mmol n'            (Opt.wwCfrac * Opt.c2n /Opt.cmw) * 1e9
    'mol N'             (Opt.wwCfrac * Opt.c2n /Opt.cmw) * 1e6
    };

area = {...
    'km^2'              1
    'm^2'               1e6
    'km^3'              Opt.depth
    'm^3'               Opt.depth * 1e9
    };
    
time = {...
    'year'              1
    'yr'                1
    'day'               365
    'd'                 365
    'sec'               86400 * 365
    's'                 86400 * 365
    };

% Parse and check unit strings

tmp = textscan(oldunit, '%s%s%s', 'delimiter', '/', 'collectoutput', true);
[m1,a1,t1] = deal(tmp{1}{:});
tmp = textscan(newunit, '%s%s%s', 'delimiter', '/', 'collectoutput', true);
[m2,a2,t2] = deal(tmp{1}{:});

m1 = validatestring(m1, mass(:,1));
m2 = validatestring(m2, mass(:,1));

a1 = validatestring(a1, area(:,1));
a2 = validatestring(a2, area(:,1));

t1 = validatestring(t1, time(:,1));
t2 = validatestring(t2, time(:,1));

% Old to new conversion factors

[~, midx] = ismember({m1 m2}, mass(:,1));
m2perm1 = mass{midx(2),2}/mass{midx(1),2};

[~, aidx] = ismember({a1 a2}, area(:,1));
a2pera1 = area{aidx(2),2}/area{aidx(1),2};

[~, tidx] = ismember({t1 t2}, time(:,1));
t2pert1 = time{tidx(2),2}/time{tidx(1),2};
 
% Field conversion factors

fac = {...
    '1/T'   1/t2pert1
    'M/A'   m2perm1/a2pera1
    'M/A/T' m2perm1/a2pera1/t2pert1
    };

%------------------------------
% Convert all tables
%------------------------------

% Group data

grp1 = {'pb', 'qb', 'emigRate', 'baRate'};  % 1/T
grp2 = {'b', 'bh'};                         % M/A
grp3 = {'dtImp', 'immig', 'emig', 'ba'};    % M/A/T

grpfac = ones(1, size(A.groupdata, 2));
grpfac(ismember(A.groupdata.Properties.VariableNames, grp1)) = fac{1,2};
grpfac(ismember(A.groupdata.Properties.VariableNames, grp2)) = fac{2,2};
grpfac(ismember(A.groupdata.Properties.VariableNames, grp3)) = fac{3,2};

A.groupdata(:,:) = num2cell(bsxfun(@times, table2array(A.groupdata), grpfac));

% Fisheries catch

A.landing(:,:) = num2cell(table2array(A.landing) .* fac{3,2}); % M/A/T
A.discard(:,:) = num2cell(table2array(A.discard) .* fac{3,2}); % M/A/T

% Stanza data

stz1 = {'BABsplit'};
stz2 = {'Btot'};

if height(A.stanzadata) > 0 && tidx(2) > 2
    warning('Stanza calculation require annual rates; use these results with caution');
end

stzfac = ones(1, size(A.stanzadata,2));
stzfac(ismember(A.stanzadata.Properties.VariableNames, grp1)) = fac{1,2};
stzfac(ismember(A.stanzadata.Properties.VariableNames, grp2)) = fac{2,2};

A.stanzadata(:,:) = num2cell(bsxfun(@times, table2array(A.stanzadata), stzfac));






