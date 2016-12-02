function EM = addpedigree(EM, prop, rowname, colname, value)
%ADDPEDIGREE Add entries to the pedigree table
%
% B = addpedigree(EM, prop, rowname, colname, value)
%
% This function provides a quick way to add values to the pedigree table
% property of an ecopathmodel object. 
%
% Input variables:
%
%   EM:         ecopathmodel object
%
%   prop:       string, property value for pedigree table, should
%               correspond to one of the table properties in A
%               ('groupdata', 'dc', 'landing', 'discard', 'df',
%               'discardFate', or 'stanzadata')
%
%   rowname:    string or cell array of strings, table row names to which
%               pedigree will be added.  If empty, it will be applied to
%               all rows of the table.
%
%   colname:    string or cell array of strings, table column (i.e.
%               variable) names to which pedigree will be added.  If empty,
%               it will be applied to all columns of the table.   
%
%   value:      pedigree value

% Copyright 2016 Kelly Kearney

% Check input

props = properties(EM);
istbl = cellfun(@(x) isa(EM.(x), 'table'), props);
tblnames = props(istbl);

if ~ischar(prop) || ~ismember(prop, tblnames)
    error('prop must correspond to an ecopathmodel table property');
end

isstr = @(x) isempty(x) || ischar(x) || (iscell(x) && all(cellfun(@ischar, x)));
if ~isstr(rowname)
    error('rowname must be a string or cell array of strings');
end
if ~isstr(colname)
    error('colname must be a string or cell array of strings');
end

if isempty(rowname)
    rowname = EM.(prop).Properties.RowNames;
end
if isempty(colname)
    colname = EM.(prop).Properties.VariableNames;
end

validateattributes(value, {'numeric'}, {'nonnegative'}, '', 'value');

% Get indices corresponding to rows and columns

[tfr, ridx] = ismember(rowname, EM.(prop).Properties.RowNames);
[tfc, cidx] = ismember(colname, EM.(prop).Properties.VariableNames);

if ~all(tfr)
    str = sprintf('%s, ', rowname{~tfr});
    error('Rows not found in table: %s', str(1:end-2));
end
if ~all(tfc)
    str = sprintf('%s, ', colname{~tfc});
    error('Columns not found in table: %s', str(1:end-2));
end

[ridx, cidx] = ndgrid(ridx, cidx);
n = numel(ridx);

% Create new table and concatenate to pedigree

New = table(repmat({prop}, n, 1), ridx(:), cidx(:), ones(n,1)*value, ...
    'VariableNames', EM.pedigree.Properties.VariableNames);

EM.pedigree = [EM.pedigree; New];