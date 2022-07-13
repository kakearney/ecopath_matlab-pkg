function [xagg, yagg] = aggregatehist(xedge, x, y, varargin)
%AGGREGATEHIST Aggregate data based on histc-like bins
%
% [xagg, yagg] = aggregatehist(xedge, x, y)
% [xagg, yagg] = aggregatehist(xedge, x, y, p1, v1, ...)
%
% Input variables:
%
%   xedge:          nbin x 1 vector (see histc) or ndim x 1 cell array of
%                   nbin# x 1 vectors (similar to hist3), of edge values
%                   for binning
%
%   x:              n x ndim array, grouping variable
%
%   y:              n x m array, another variable to be grouped along with
%                   x 
%
% Optional input variables (pass as parameter/value pairs)
%
%   combineupper:   logical scalar, if true, combine the upper-edge bin
%                   with the previous bin, such that there are
%                   nbin-1 bins.  If false, follows same rules as
%                   in histc, resulting in nbin bins.
%
% Output variables:
%
%   xagg:           nbin1 x nbin2 x ... or nbin1-1 x nbin2-1 x ... cell
%                   array of aggregated x values 
%
%   yagg:           nbin1 x nbin2 x ... or nbin1-1 x nbin2-1 x ... cell
%                   array of aggregated y values 

% Copyright 2012 Kelly Kearney

% Check input

if ndims(y) > 2
    error('Not sure this works for higher-dimensional arrays');
end

if isvector(x)
    x = x(:);
end

[nrow, ncol] = size(x);
if ncol == 1 && ~iscell(xedge)
    xedge = {xedge};
end

if size(y,1) ~= nrow
    error('x and y must have same number of rows');
end

p = inputParser;
p.addParameter('combineupper', true, @(x) validateattributes(x,{'logical'}, {'scalar'}));
p.parse(varargin{:});
Opt = p.Results;

% Opt.combineupper = true;
% Opt = parsepv(Opt, varargin);

% Bin x data

[bin, sz] = ndhistc(x, xedge);
sz = reshape(sz, 1, []); % Needs to be row vector

if any(bin == 0)
    warning('aggregatehist:outside', 'Some data outside bin edges or NaN');
end

% Aggregate x and y data

[idx, xtmp] = aggregate(bin, x);
[idx, ytmp] = aggregate(bin, y);

xagg = cell(sz);
yagg = cell(sz);

isin = idx ~= 0;

xagg(idx(isin)) = xtmp(isin);
yagg(idx(isin)) = ytmp(isin);
        
% Combine upper bins if specifies

if Opt.combineupper
    if size(x,2) == 1
        xcombo = [xagg{end-1}; xagg{end}];
        xagg{end-1} = xcombo;
        xagg = xagg(1:end-1);
        
        ycombo = [yagg{end-1}; yagg{end}];
        yagg{end-1} = ycombo;
        yagg = yagg(1:end-1);
        
    elseif size(x,2) == 2
        xcombo = cellfun(@(a,b) [a;b], xagg(end-1,:), xagg(end,:), 'uni', 0);
        xagg(end-1,:) = xcombo;
        xagg = xagg(1:end-1,:);
        xcombo = cellfun(@(a,b) [a;b], xagg(:,end-1), xagg(:,end), 'uni', 0);
        xagg(:,end-1) = xcombo;
        xagg = xagg(:,1:end-1);
        
        ycombo = cellfun(@(a,b) [a;b], yagg(end-1,:), yagg(end,:), 'uni', 0);
        yagg(end-1,:) = ycombo;
        yagg = yagg(1:end-1,:);
        ycombo = cellfun(@(a,b) [a;b], yagg(:,end-1), yagg(:,end), 'uni', 0);
        yagg(:,end-1) = ycombo;
        yagg = yagg(:,1:end-1);
        
    else
        error('Combineupper option not supported for ndims>3');
    end
end


% Subfunction: n-dimensional histogram, with bin #s returned

function [binidx, sz] = ndhistc(xval, xedge)

sz = cellfun(@length, xedge);
if isscalar(sz)
    sz = [sz 1];
end

bin = zeros(size(xval));
for ix = 1:size(xval,2)
    [blah, bin(:,ix)] = histc(xval(:,ix), xedge{ix});
end
isout = any(bin==0, 2);
bin = num2cell(bin(~isout,:), 1);
binidx = zeros(size(xval,1), 1);
binidx(~isout) = sub2ind(sz, bin{:});




