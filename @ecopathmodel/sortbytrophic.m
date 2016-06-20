function [A, isrt] = sortbytrophic(A)
%SORTBYTROPHIC Sort ecopathmodel object groups by trophic level
%
% [A, isrt] = sortbytrophic(A)
%
% Input variables:
%
%   A:      ecopathmodel object
%
% Output variables:
%
%   B:      ecopathmodel object, same as A but with groups arranged in
%           descending order by trophic level
%
%   isrt:   indices corresponding to location of each output model group in
%           the input model 

% First, check to make sure detritus groups come last.  If not, the trophic
% level calculations will be wrong

names = A.name;

[~, isrtpp] = sort(A.groupdata.pp);
A = A.sort(isrtpp);

% Now sort by trophic level

tl = trophiclevel(A);
[~, isrttl] = sort(tl, 'descend');

A = A.sort(isrttl);

% Indices

[~,isrt] = ismember(A.name, names);