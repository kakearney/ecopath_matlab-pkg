function B = subpedigreevalues(A, vals)
%SUBPEDIGREEVALUES Replace values in ecopathmodel based on pedigree
%
% B = subpedigreevalues(A, vals)
%
% This function returns matrices, corresponding to the table in the
% ecopathmodel object A, with new values subtituted where indicated by the
% pedigree table.
%
% Input variables:
%
%   A:      ecopathmodel object.
%
%   vals:   nped x n array of values, where rows of this matrix correspond
%           to the rows in the A.pedigree table, and columns each represent
%           a single set of values to substitute into A
%
% Output variables:
%
%   B:      structure with fields corresponding to tables in A.  Each field
%           holds a nrow x ncol x n array, where nrow and ncol correspond
%           to the numbers of rows and columns in the corresponding table
%           in A, and n is the number of new parameter sets.  Parameter
%           sets will be adjusted to meet constraints (for example, diet
%           components will be renormalized so predator diet sums to 1, and
%           parameters related to the stable growth curve of stanza groups
%           must be synced up).

% Copyright 2016 Kelly Kearney

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
    B.dc = bsxfun(@rdivide, B.dc, sum(B.dc,1));
    B.dc(isnan(B.dc)) = 0;
end
    
% If any of the parameters involved in the stable age growth curve were
% changed, recalculate the values that depend on this curve.

idx = stanzaindices(A);

[~,loc] = ismember({'ageStart', 'vbK', 'b', 'qb', 'pb', 'ba'}, ...
    A.groupdata.Properties.VariableNames);

for is = 1:length(A.stanza)
    usebtot = isnan(A.groupdata.b(idx{is}(end))) & ~isnan(A.stanzadata.Btot(is));
    
    if isfield(B, 'groupdata')
        a = permute(B.groupdata(idx{is},   loc(1),:), [1 3 2]);
        k = permute(B.groupdata(idx{is}(1),loc(2),:), [1 3 2]);
        if ~usebtot
            blead = permute(B.groupdata(idx{is}(end), loc(3),:), [1 3 2]);
        end
        qblead = permute(B.groupdata(idx{is}(end), loc(4),:), [1 3 2]);
        z = permute(B.groupdata(idx{is}, loc(5),:), [1 3 2]);
    else
        a = repmat(A.groupdata.ageStart(idx{is}), 1, nnew);
        k = repmat(A.groupdata.vbK(idx{is}(1)), 1, nnew);
        if ~usebtot
            blead = repmat(A.groupdata.b(idx{is}(end)), 1, nnew);
        end
        qblead = repmat(A.groupdata.qb(idx{is}(end)), 1, nnew);
        z = repmat(A.groupdata.pb(idx{is}), 1, nnew);
    end
    
    [~,sloc] = ismember({'BaBsplit', 'Btot'},  A.stanzadata.Properties.VariableNames);
    if isfield(B, 'stanzadata')
        bab = permute(B.stanzadata(is,sloc(1),:), [1 3 2]);
        if usebtot
            blead = permute(B.stanzadata(is,sloc(2),:), [1 3 2]);
        end
    else
        bab = repmat(A.stanzadata.BABsplit(is), 1, nnew);
        if usebtot
            blead = repmat(A.stanzadata.Btot(is), 1, nnew);
        end
    end
    
    allvars = [a; k; blead; qblead; z; bab];
    changed = any(bsxfun(@minus, allvars, allvars(:,1)), 2) & ...
             ~all(isnan(allvars),2);
    
    if usebtot
        bflag = 'total';
    else
        bflag = 'lead';
    end
         
    if any(changed)
        if ~isfield(B, 'groupdata')
            B.groupdata = repmat(table2array(A.groupdata), 1, 1, nnew);
        end
        for in = 1:nnew
            Out = editstanzacalcs(a(:,in), k(:,in), bab(:,in), ...
                blead(:,in), qblead(:,in), z(:,in), 1, bflag); 
            B.groupdata(idx{is},loc(3),in) = Out.b;
            B.groupdata(idx{is},loc(4),in) = Out.qb;
            B.groupdata(idx{is},loc(6),in) = Out.ba;
        end
    end
    
end






