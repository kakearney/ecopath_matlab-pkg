function New = combinegroups(EM, varargin)
%COMBINEGROUPS Combine groups and/or fleets in an ecopathmodel object
%
% New = combinegroups(EM, grp1, grp2, ...)
% New = combinegroups(EM, grp1, grp2, 'labels', labels)
%
% This function consolidates two or more functional groups (or fleets) into
% one group, attempting to preserve the magnitude of fluxes into and out of
% the groups.  
%
% Groups to be combined must be of the same type (for example, you cannot
% combine a producer with a consumer or a detrital group).  
%
% Input variables:
%
%   EM:     ecopathmodel object
%
%   grp#:   cell array of strings.  Each cell array represents a list of
%           functional groups or fleets that are to be combined.
%
%   labels: n x 1 cell array of strings, where n is the number of grp cell
%           arrays passed as input.  These labels are the new group/fleet
%           names used for the combined groups.  If not included, new
%           groups will be named newgroup1, newgroup2, etc.
%
% Output variables:
%
%   New:    ecopathmodel object with new consolidated groups

% Copyright 2016 Kelly Kearney

%--------------------------
% Setup
%--------------------------

% Check input
% TODO: make sure not combining mixed-type, mixed-stanza-set,
% non-contiguous stanza groups, etc

isname = cellfun(@(x) ischar(x) && strcmp(x,'labels'), varargin);
if any(isname)
    idx = find(isname);
    labels = varargin{idx+1};
    isname(idx+1) = true;
    grps = reshape(varargin(~isname), [], 1);
else
    grps = varargin(:);
    labels = cellstr(num2str((1:length(grps))', 'newgroup%d'));
end
labels = labels(:);
grps = grps(:);

iscellstr = @(x) iscell(x) && all(cellfun(@ischar, x));

if ~iscell(grps) || ~all(cellfun(iscellstr, grps))
    error('Inputs should be cell arrays of strings');
end

if ~(iscellstr(labels) && isequal(size(grps), size(labels)))
    error('labels array should be cell array of strings, with same # of elements as groups listed');
end

grps = cellfun(@(x) x(:), grps, 'uni', 0); % make sure all column arrays

tf = ismember(cat(1, grps{:}), [EM.name; EM.fleet]);
if ~all(tf)
    tmp = cat(1, grps{:});
    str = sprintf('%s, ', tmp{~tf});
    error('Input name(s) doesn''t match model group/fleet names: %s', str(1:end-2));
end

% Aggregating indices:
% Ordered by grouped stuff, then remaining groups, then fleets

% Calculate new-group index array

gidx = (1:length(grps))';

idx = cellfun(@(x,n) ones(size(x))*n, grps, num2cell(gidx), 'uni', 0);
idx = cat(1, idx{:});

[tf, loc] = ismember([EM.name; EM.fleet], cat(1, grps{:}));

aidx = nan(EM.ngroup + EM.ngear,1);
aidx(tf) = idx(loc(tf));

aidx(~tf) = (1:sum(~tf)) + max(gidx);

oldnames = [EM.name; EM.fleet];
oldnames = [labels(:); oldnames(~tf)];
newnames = oldnames(aidx);

[~,~,aidx] = unique(aidx, 'stable'); % Renumber to keep live, det, and fleets in same order

aidxg = aidx(1:EM.ngroup);
aidxd = aidx((EM.nlive+1):EM.ngroup);
aidxf = aidx((1:EM.ngear)+EM.ngroup);

impidx = max(aidx) + 1; % import
expidx = impidx + 1;    % export

% Names of new functional groups and fleets

[~, newgrpname] = aggregate(aidxg, newnames(1:EM.ngroup), @(x) x{1});
[~, newfltname] = aggregate(aidxf, newnames(EM.ngroup+(1:EM.ngear)), @(x) x{1});

% Types and stanzas for each combined set (check that not mixing types or
% stanzas)

[~, pptype] = aggregate(aidx, [EM.groupdata.pp; ones(EM.ngear,1)*3]);
if any(cellfun(@(x) length(unique(x)), pptype) > 1)
    error('Cannot combine groups of different types');
end
pptype = cellfun(@(x) x(1), pptype);

[~, stz] = aggregate(aidxg, EM.groupdata.stanza);
if any(cellfun(@(x) length(unique(x)), stz) > 1)
    error('Cannot combine groups from different stanza sets');
end

% Check that only contiguous-age stanza groups are being combined

for ii = 1:length(stz)
    if stz{ii}(1) > 0 % includes multi-stanza groups
       grpsinset = find(EM.groupdata.stanza == stz{ii}(1));
       grpsincombo = find(aidxg == ii);
       [~,ageorder] = sort(EM.groupdata.ageStart(grpsinset));
       grpsinset = grpsinset(ageorder);
       [~,loc] = ismember(grpsincombo,grpsinset);
       if ~all(diff(sort(loc)) == 1)
           error('Only continguos-age stanza groups can be combined');
       end
    end
end

% Sum biomass (use balanced version)

[Ep, epflag] = EM.ecopath;

[~, b] = aggregate(aidxg, Ep.b);

bsum = cellfun(@sum, b);
bfrac = cellfun(@(x,y) x./y, b, num2cell(bsum), 'uni', 0);

%--------------------------
% Aggregate properties
%--------------------------

% Diet

dietflux = [Ep.q0; Ep.import(1:EM.ngroup)'];
dietflux(:,(EM.nlive+1):EM.ngroup) = 0; % remove flow to detritus and detritus import

dietagg = aggregate2d(dietflux, [aidxg; impidx], aidxg);
dc = bsxfun(@rdivide, dietagg, sum(dietagg,1));
dc(isnan(dc)) = 0;

import = dc(end,:);
dc = dc(1:end-1,:);

% ageStart, minimum of groups

[~, age] = aggregate(aidxg, EM.groupdata.ageStart, @min);
age = cat(1, age{:});

% Detritus import, sum

[~, dtimp] = aggregate(aidxg, EM.groupdata.dtImp, @sum);
dtimp = cat(1, dtimp{:});

% Immigration/emigration

[~,imig] = aggregate(aidxg, EM.groupdata.immig, @sum);
[~,emig] = aggregate(aidxg, EM.groupdata.emig, @sum);
imig = cat(1, imig{:});
emig = cat(1, emig{:});

% Main Ecopath variables: back-calculate from Ecopath balance

P = Ep.b.*Ep.pb;
Q = Ep.qb.*Ep.b;
M0 = Ep.b.*Ep.otherMortRate;
R = Ep.respiration;

[~, P] = aggregate(aidxg, Ep.b.*Ep.pb, @sum);
[~, Q] = aggregate(aidxg, Ep.b.*Ep.qb, @sum);
[~,M0] = aggregate(aidxg, Ep.b.*Ep.otherMortRate, @sum);
[~, R] = aggregate(aidxg, Ep.respiration, @sum);
[~, Y] = aggregate(aidxg, sum(table2array(EM.landing) + table2array(EM.discard),2), @sum);
[~,M2] = aggregate(aidxg, Ep.predMortRate.*Ep.b, @sum);

P  = cell2mat(P);
Q  = cell2mat(Q);
M0 = cell2mat(M0);
R  = cell2mat(R);
Y  = cell2mat(Y);
M2 = cell2mat(M2);

E = emig - imig;

idxout = max(aidx) + (1:4)';
flows = aggregate2d(Ep.flow, [aidx; idxout], [aidx; idxout]);
Idx.det = find(pptype == 2);
Idx.res = size(flows,1) - 2;
Idx.liv = find(pptype < 2);
detin = flows(:,Idx.det);
deteaten = flows(Idx.det,Idx.liv);
detresp = flows(Idx.det, Idx.res);
detee = sum(deteaten,2)'./(sum(detin,1) - sum(detresp,2)');

ngnew = length(bsum);

gd.b = bsum;
gd.pb = P./bsum;
gd.qb = Q./bsum;
gd.ee = 1 - M0./P; % except det
gd.ee(pptype == 2) = detee;
gd.ge = P./Q;
gd.ge(pptype(1:ngnew)>=1) = 0;
gd.gs = (Q - P - R)./Q;
gd.gs(pptype(1:ngnew)>=1) = 0;

[~,hasba] = aggregate(aidxg, EM.groupdata.ba, @(x) any(x~=0));
hasba = cell2mat(hasba);
gd.ba = zeros(size(gd.b));
if any(hasba) % b/c floating point, only do if necessary
    batmp = P - Y - M2 - E - P.*(1-gd.ee);
    gd.ba(hasba) = batmp(hasba);    
end

[~, detpb] = aggregate(aidxg, EM.groupdata.detpb, @nanmean);
detpb = cat(1, detpb{:});

% Add back the specially-handled groupdata variables

gd.ageStart = age;
gd.import = import';
gd.dtImp = dtimp;
gd.b = bsum;
gd.pp = pptype(1:ngnew);
gd.areafrac = ones(ngnew,1);
gd.immig = imig;
gd.emig = emig;
gd.detpb = detpb;

% Rate inputs

gd.baRate = gd.ba./gd.b;
gd.emigRate = gd.emig./gd.b;


% If parameter was missing from all sub-groups in the original model,
% remove it again (it there was a mix, keep it)

gdvar = EM.groupdata.Properties.VariableNames;

tmp = nan(length(gd.b), length(gdvar));
for ii = 1:length(gdvar)
    if isfield(gd, gdvar{ii})
        tmp(:,ii) = gd.(gdvar{ii});
    end
end

[~, gdmask] = aggregate(aidxg, ismissing(EM.groupdata), @(x) all(x,1));
gdmask = cat(1, gdmask{:});

tmp(gdmask) = NaN;
gd = array2table(tmp, 'variablenames', EM.groupdata.Properties.VariableNames);

% Landing and discard: sum

landing = aggregate2d(table2array(EM.landing), aidxg, aidxf);
discard = aggregate2d(table2array(EM.discard), aidxg, aidxf);

% Detritus and discard fate

detf = [Ep.q0(:,(EM.nlive+1):EM.ngroup) Ep.export(1:EM.ngroup)]; 
disf = [Ep.flow(Ep.Idx.gear, Ep.Idx.det) Ep.export(EM.ngroup+(1:EM.ngear))];

detfagg = aggregate2d(detf, aidxg, [aidxd; expidx]);
disfagg = aggregate2d(disf, aidxf, [aidxd; expidx]);

detf = bsxfun(@rdivide, detfagg, sum(detfagg,2));
disf = bsxfun(@rdivide, disfagg, sum(disfagg,2));
disf(isnan(disf)) = 0;  % TODO: if group has no discard, this doesn't preserve discard fate... but does that matter?

detf = detf(:,1:end-1);
disf = disf(:,1:end-1);


% Stanzas

stz = cellfun(@(x) x(1), stz);
keepstz = false(height(EM.stanzadata),1);
for is = 1:height(EM.stanzadata)
    ns = sum(stz == is);
    keepstz(is) = ns > 1;
    if ~keepstz(is)
        stz(stz == is) = 0;
    end
end

sidxold = 1:height(EM.stanzadata);
sidxold = sidxold(keepstz);
[~,~,sidxnew] = unique(sidxold);
[tf,loc] = ismember(stz, sidxold);
stz(tf) = sidxnew(loc(tf));

sd = EM.stanzadata(keepstz,:);
sd.stanzaID = sidxnew;

[~,vbk] = aggregate(aidxg, EM.groupdata.vbK, @(x) x(1));
vbk = cat(1, vbk{:});

gd.stanza = stz;
gd.vbK(stz > 0) = vbk(stz > 0);

% Pedigree: straight average

hasped = ~isempty(EM.pedigree);

if hasped
    [tbl, ped] = aggregate(EM.pedigree.property, EM.pedigree);
    pedprop = cell(0);
    pedval = zeros(0,3);
    for ii = 1:length(tbl)
        switch tbl{ii}
            case 'dc'
                rnew = aidxg(ped{ii}.row);
                cnew = aidxg(ped{ii}.column);
            case 'discard'
                rnew = aidxg(ped{ii}.row);
                [~,~,funq] = unique(aidxf, 'stable');
                cnew = funq(ped{ii}.column);
            case 'groupdata'
                rnew = aidxg(ped{ii}.row);
                cnew = ped{ii}.column;
            case 'landing'
                rnew = aidxg(ped{ii}.row);
                [~,~,funq] = unique(aidxf, 'stable');
                cnew = funq(ped{ii}.column);
            otherwise
                error('Combining pedigree for properties other than groupdata, dc, discard, or landing not yet supported'); 
        end

        if isscalar(rnew)
            rc = [rnew cnew];
            pd = ped{ii}.pedigree;
        else
            [rc, pd] = aggregate([rnew cnew], ped{ii}.pedigree, @mean);
            pd = cat(1, pd{:});
        end
            

        nnew = size(rc,1);
        pedprop = [pedprop; repmat(tbl(ii), nnew, 1)];
        pedval = [pedval; [rc pd]];

    end

    Ped = table(pedprop, pedval(:,1), pedval(:,2), pedval(:,3), 'VariableNames', EM.pedigree.Properties.VariableNames);
end
%--------------------------
% New object
%--------------------------

New = ecopathmodel(ngnew, sum(pptype < 2), sum(pptype == 3), ...
    'pp', pptype(1:ngnew), 'groups', newgrpname, 'fleets', newfltname);

New.stanzadata = sd;
New.stanza = sd.Properties.RowNames;

for ii = 1:length(gdvar)
    New.groupdata.(gdvar{ii}) = gd.(gdvar{ii});
end

New.dc(:,:) = num2cell(dc);
New.landing(:,:) = num2cell(landing);
New.discard(:,:) = num2cell(discard);
New.df(:,:) = num2cell(detf);
New.discardFate(:,:) = num2cell(disf);

if hasped
    New.pedigree = Ped;
end

New = New.sortbytrophic;





% Subfunctions ************************************************************

    

function xgrp = aggregate2d(x, ridx, cidx)

[~, xtmp] = aggregate(ridx, x, @(x) sum(x,1));
xtmp = cat(1, xtmp{:});

[~, xtmp] = aggregate(cidx, xtmp', @(x) sum(x,1));
xgrp = cat(1, xtmp{:})';

