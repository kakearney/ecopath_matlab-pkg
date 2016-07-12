function combinegroups(EM, grps)

% Calculate new-group index array

grps = grps(:);
gidx = (1:length(grps))';

idx = cellfun(@(x,n) ones(size(x))*n, grps, num2cell(gidx), 'uni', 0);
idx = cat(1, idx{:});

[tf, loc] = ismember(EM.name, cat(1, grps{:}));

aidx = nan(EM.ngroup,1);
aidx(tf) = idx(loc(tf));

aidx(~tf) = (1:sum(~tf)) + max(gidx);

% Aggregate biomass (use balanced version)

Ep = EM.ecopath;

[~, b] = aggregate(aidx, Ep.b);

bsum = cellfun(@sum, b);
bfrac = cellfun(@(x,y) x./y, b, num2cell(bsum), 'uni', 0);

% Calculate weighted-average values for all other parameters



blah