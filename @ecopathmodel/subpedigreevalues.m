function B = subpedigreevalues(A, vals)
%SUBPEDIGREEVALUES Replace values in ecopathmodel based on pedigree


validateattributes(A, {'ecopathmodel'}, {'scalar'}, 1);
nped = height(A.pedigree);
validateattributes(vals, {'numeric'}, {'2d', 'nrows', nped});

nnew = size(vals,2);

% First, simply substitute the supplied values into the appropriate tables

[tbl, ijval] = aggregate(A.pedigree.property, ...
    [A.pedigree.row A.pedigree.column vals]);

B = struct;

for it = 1:length(tbl)
    B.(tbl{it}) = repmat(table2array(A.(tbl{it})), 1, 1, nnew);
    
    ridx = repmat(ijval{it}(:,1), 1, nnew);
    cidx = repmat(ijval{it}(:,2), 1, nnew);
    pidx = repmat(1:nnew, size(ijval{it},1), 1);
    idx = sub2ind(size(B.(tbl{it})), ridx, cidx, pidx);
    B.(tbl{it})(idx) = ijval{it}(:,3:end);
end

% If the diet table was altered, make sure the values sum to 1

if isfield(B, 'dc')
    B.dc = bsxfun(@rdivide, B.dc, sum(B.dc,2));
    B.dc(isnan(B.dc)) = 0;
end
    
% If the B, PB (i.e. Z), QB, or BA of a leading stanza group was
% altered, or the Btot or BaB value for a stanza set, recalculate
% stable age distribution. 


