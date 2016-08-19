function [xcon, yagg, yidxagg] = aggregate(x, y, fun)
%AGGREGATE Aggregate values into cell array
%
% [xcon, yagg, yidxagg] = aggregate(x, y)
% [xcon, yagg, yidxagg] = aggregate(x, y, fun)
%
% This function groups together values of y, based on category values in x.
% It performs more or less like accumaray(x,y,[a b]. @(x) {x}), except
% allows x to be any value, not just indices, and y can have any number of
% columns.  
%
% Input variables:
%
%   x:          n x 1 array, categories, can be either numeric or a cell
%               array of strings  
%
%   y:          n x m array, values to be grouped
%
%   fun:        function handle. If included, this function is applied to
%               the grouped values of y 
%
% Output variables:
%
%   xcon:       unique values of x
%
%   yagg:       cell array of y values corresponding to each x.
%
%   yidxagg:    row indices of aggregated values

% Copyright 2013 Kelly Kearney

if ~isequal(size(x,1), size(y,1))
    error('x and y must have same number of columns');
end

if nargin > 2
    applyfun = true;
else
    applyfun = false;
end

% Get x values accumarray will take

if isvector(x)
    [xcon, blah, ix] = unique(x);
elseif isnumeric(x)
    [xcon, blah, ix] = unique(x, 'rows');
elseif iscell(x)
    loc = zeros(size(x));
    xunq = cell(size(x,2),1);
    for ic = 1:size(x,2)
        xunq{ic} = unique(x(:,ic));
        [tf, loc(:,ic)] = ismember(x(:,ic), xunq{ic});
    end
    [loccon, blah, ix] = unique(loc, 'rows');
    xcon = cell(size(loccon));
    for ic = 1:size(x,2)
        xcon(:,ic) = xunq{ic}(loccon(:,ic));
    end
end

nrow = max(ix);
ncol = size(y,2);

% Accumarray the y indices and translate that to y

yrow = (1:size(y,1))';
yidxagg = accumarray(ix, yrow, [nrow 1], @(x) {x});

yagg = cell(size(yidxagg));
for iy = 1:length(yagg)
    if applyfun
        yagg{iy} = fun(y(yidxagg{iy},:));
    else
        yagg{iy} = y(yidxagg{iy},:);
    end 
end



