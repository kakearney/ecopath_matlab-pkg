function [EM,B] = rpath2ecopathmodel(basename, varargin)
%RPATH2ECOPATHMODEL Create ecopathmodel object from Rpath data files
%
% [A, B] = rpath2ecopathmodel(basename, p1, v1, ...)
%
% Rpath (the R-based version of Ecopath/Ecosim, available from
% https://github.com/slucey/RpathDev) allows model data be imported from
% and exported to .csv files.  This function reads data from a set of these
% files into an ecopathmodel object.
% 
% Input variables:
%
%   basename:   Base file name for the 4 input files
%
% Optional input variables (passed as paramete/value pairs)
%
%   basestr:    string appended to basename for basic data, not including
%               .csv extension ['_model']
%
%   dietstr:    string appended to basename for diet data, not including
%               .csv extension ['_diet']
%
%   juvsstr:    string appended to basename for old-style stanza data (only 
%               applicable if old = true), not including .csv extension
%               ['_juvs'] 
%
%   pedstr:     string appended to basename for pedigree data, not
%               including .csv extension ['_pedigree'] 
%
%   stanzastr:  string appended to basename for stanza data (only 
%               applicable if old = false), not including .csv extension
%               ['_stanzas']
%
%   stgrpstr:   string appended to basename for stanza group data (only 
%               applicable if old = false), not including .csv extension
%               ['_staanza_groups']
%
%   old:        logical scalar, true if model has Kerim's old single-table
%               _juvs file rather than the newer 2-table stanza and
%               stanza_groups files. [false]
%
% Output variables:
%
%   EM:         ecopathmodel object
%
%   B:          structure with full input tables from files

%--------------------
% Parse input
%--------------------

p = inputParser;
p.addParameter('basestr', '_model',    @(x) validateattributes(x, {'char'}, {}));
p.addParameter('dietstr', '_diet',     @(x) validateattributes(x, {'char'}, {}));
p.addParameter('juvsstr', '_juvs',     @(x) validateattributes(x, {'char'}, {}));
p.addParameter('pedstr',  '_pedigree', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('stanzastr', '_stanzas', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('stgrpstr', '_stanza_groups', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('old', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.parse(varargin{:});

Opt = p.Results;


%--------------------
% Read tables
%--------------------

% Read file data into tables

Base = readtable([basename Opt.basestr '.csv'], 'TreatAsEmpty', 'NA');
Diet = readtable([basename Opt.dietstr '.csv'], 'ReadRowNames', true, 'TreatAsEmpty', 'NA');
Ped  = readtable([basename Opt.pedstr '.csv'], 'TreatAsEmpty', 'NA');

if Opt.old
    Juvs = readtable([basename Opt.juvsstr '.csv']);
   
    ns = height(Juvs);
    
    Sgrp = table((1:ns)', Juvs.StanzaName, ones(ns,1)*2, Juvs.VonBK, ...
        Juvs.VonBD, Juvs.Wmat50, Juvs.RecPower, 'VariableNames', ...
        {'StGroupNum', 'StanzaGroup', 'nstanzas', 'VBGF_Ksp', 'VBGF_d', ...
        'Wmat', 'RecPower'});
    
    grpnum = [Juvs.JuvNum Juvs.AduNum]';
    grpnum = grpnum(:);
    
    nsg = length(grpnum);
    
    z = Base.PB(grpnum) - Base.BioAcc(grpnum); % Note: this is prob. wrong, but it's just a placeholder for now (I don't use it; Rpath does).
    
    Stan = table(kron((1:ns)', ones(2,1)), repmat([1;2],ns,1), grpnum, ...
        Base.Group(grpnum), nan(nsg,1), nan(nsg,1), z, false(nsg,1), ...
        'VariableNames', {'StGroupNum', 'Stanza', 'GroupNum', 'Group', ...
        'First', 'Last', 'Z', 'Leading'});
    Stan.First(Stan.Stanza == 1) = 0;
    Stan.First(Stan.Stanza == 2) = Juvs.RecAge*12;
    
    Stan.Last(Stan.Stanza == 1) = Juvs.RecAge*12 - 1;
    Stan.Last(Stan.Stanza == 2) = 999;
    
    Stan.Leading = Stan.Stanza == 2;
    
else
    Stan = readtable([basename Opt.stanzastr '.csv'], 'TreatAsEmpty', 'NA');
    Sgrp = readtable([basename Opt.stgrpstr '.csv'], 'TreatAsEmpty', 'NA');
 
    isn = isnan(Stan.StGroupNum);
    Stan = Stan(~isn,:);
    
    isn = isnan(Sgrp.StGroupNum);
    Sgrp = Sgrp(~isn,:);
    
end

B.Base = Base;
B.Diet = Diet;
B.Ped = Ped;
B.Stan = Stan;
B.Sgrp = Sgrp;
if Opt.old
    B.Juvs = Juvs;
else
    B.Juvs = [];
end

sflag = ~isempty(Stan);

% Parse column names

bcol = parsecolname(Base);
dcol = parsecolname(Diet);
pcol = parsecolname(Ped);

isgroup = Base.Type <= 2;
isgear  = Base.Type == 3;
isdet   = Base.Type == 2;

%--------------------
% Build ecopathmodel
%--------------------

ngroup = sum(isgroup);
nlive  = ngroup - sum(isdet);
ngear  = sum(isgear);
pp     = Base.Type(isgroup);

groups = strtrim(Base.Group(isgroup));
fleets = strtrim(Base.Group(isgear));
if sflag
    stanzas = Sgrp.StanzaGroup;
else
    stanzas = {};
end

grp = alternames(groups);
flt = alternames(fleets);
stz = alternames(stanzas);


old = [groups; fleets; stanzas];
new = [grp; flt; stz];
issame = strcmp(old, new);
if ~all(issame)
    tmp = [old(~issame) new(~issame)]';
    str = sprintf('  %s -> %s\n', tmp{:});
    warning('Ecopathmodel:rpath2ecopathmodel:names', 'Some names changed to meet Matlab''s variable name restrictions:\n%s', str);
end

EM = ecopathmodel(ngroup, nlive, ngear, 'groups', grp, 'fleets', flt, ...
    'stanzas', stz, 'pp', pp);

% Fill in basic data

EM.groupdata.areafrac = ones(ngroup,1);
EM.groupdata.b        = Base.Biomass(isgroup);
EM.groupdata.pb       = Base.PB(isgroup);
EM.groupdata.qb       = Base.QB(isgroup);
EM.groupdata.ee       = Base.EE(isgroup);
EM.groupdata.ge       = Base.ProdCons(isgroup);
EM.groupdata.gs       = Base.Unassim(isgroup);
EM.groupdata.dtImp    = Base.DetInput(isgroup);
EM.groupdata.ba       = Base.BioAcc(isgroup);

EM.groupdata.immig    = zeros(ngroup,1);
EM.groupdata.emig     = zeros(ngroup,1);

% Diet

[tfr, locr] = ismember(Diet.Properties.RowNames, groups);
[tfc, locc] = ismember(dcol, groups);

dc = zeros(ngroup);
dc(locr(tfr), locc(tfc)) = table2array(Diet(tfr,tfc));
dc(isnan(dc)) = 0;

EM.dc(:,:) = num2cell(dc);

% Flow to detritus

ndet = ngroup - nlive;

idxdet = (1:ndet)+10;
idxfish = (1:ngear)+max(idxdet);
idxdisc = (1:ngear)+max(idxfish);

EM.df(:,:) = Base(isgroup,idxdet);

% Fisheries

EM.landing(:,:) = Base(isgroup,idxfish);
EM.discard(:,:) = Base(isgroup,idxdisc);

EM.discardFate(:,:) = Base(isgear,idxdet);

% Stanza data

EM.groupdata.stanza(Stan.GroupNum) = Stan.StGroupNum;
EM.groupdata.ageStart(Stan.GroupNum) = Stan.First;

[tf, loc] = ismember(EM.groupdata.stanza, Sgrp.StGroupNum);

EM.groupdata.vbK(tf) = Sgrp.VBGF_Ksp(loc(tf));

EM.stanzadata.stanzaID = Sgrp.StGroupNum;

% The way Rpath deals with BA across stanza sets has changed over time, and
% still isn't completely clear to me... some of Kerim's tables include
% different values in the _juvs and _base files. Consider this bit a work
% in progress still.

if Opt.old
    % Allows different BA for juvs and adults, and doesn't always match BA
    bab1 = [Juvs.JuvZ_BAB, Juvs.AduZ_BAB];
    ajidx = [Juvs.JuvNum Juvs.AduNum];
    bab2 = Base.BioAcc(ajidx)./Base.Biomass(ajidx);
    
    if ~isequal(bab1, bab2)
        tmp = [EM.name(ajidx(:)) num2cell([bab2(:) bab1(:)])]';
        str = sprintf('  %s: %g (base) vs %g (juvs)\n', tmp{:});
        warning('Found conflicting BA/B data for multi-stanza groups in Base and Juvs tables:\n%sReplacing Base values with those from Juvs', str);
    end
    
    if isequal(bab1(:,1), bab1(:,2))
        EM.stanzadata.BABsplit = bab1(:,1);
    else
        EM.stanzadata.BABsplit(:) = NaN;
        EM.groupdata.ba([Juvs.JuvNum; Juvs.AduNum]) = NaN;
        EM.groupdata.baRate([Juvs.JuvNum; Juvs.AduNum]) = bab1(:);
    end

%     EM.stanzadata.BABsplit = arrayfun(@(x,y) [x y], Juvs.JuvZ_BAB, Juvs.AduZ_BAB, 'uni', 0);
else
    % Will have to ask Sean how this is stored now.
end

% Pedigree: Rpath allows pedigree values for B, PB, QB, DC, and each
% fishing fleet, on a per-group (or per-caught-group) basis.  For diet, the
% same value is applied to all prey components of a group's diet.  For
% catch, the same value is applied to landing and discard of each pairing.  

[~,loc] = ismember({'b','pb','qb','dc'}, EM.groupdata.Properties.VariableNames);

pedprop = repmat({'groupdata'}, ngroup*3,1);
pedrow = repmat((1:ngroup)', 3, 1);
pedcol = kron(loc(1:3)', ones(ngroup,1));
pedval = [Ped.B(1:ngroup); Ped.PB(1:ngroup); Ped.QB(1:ngroup)];

dc = table2array(EM.dc);
dcped = repmat(Ped.Diet(1:ngroup)', ngroup, 1);
[ipry, iprd] = find(dc);

pedprop = [pedprop; repmat({'dc'}, length(ipry), 1)];
pedrow = [pedrow; ipry];
pedcol = [pedcol; iprd];
pedval = [pedval; dcped(dc > 0)];

gearped = table2array(Ped(1:ngroup,5:end));
land = table2array(EM.landing);
disc = table2array(EM.discard);
[ipry,iflt] = find(land);

pedprop = [pedprop; repmat({'landing'}, length(ipry), 1)];
pedrow = [pedrow; ipry];
pedcol = [pedcol; iflt];
pedval = [pedval; gearped(land > 0)];

[ipry,iflt] = find(disc);

pedprop = [pedprop; repmat({'discard'}, length(ipry), 1)];
pedrow = [pedrow; ipry];
pedcol = [pedcol; iflt];
pedval = [pedval; gearped(disc > 0)];

EM.pedigree = table(pedprop, pedrow, pedcol, pedval, ...
    'VariableNames', {'property', 'row', 'column', 'pedigree'});



%--------------------
% Subfunctions
%--------------------

% Get original table column names, rather than the altered ones

function name = parsecolname(T)
name = T.Properties.VariableNames;
tmp = regexp(T.Properties.VariableDescriptions, 'Original column heading: ''(.*)''', 'tokens', 'once');
isemp = cellfun('isempty', tmp);
tmp = cat(1, tmp{:});
name(~isemp) = tmp;


function x = alternames(x)

isgood = cellfun(@isvarname, x);

% Often named using ages, so move numeric bits to end

x(~isgood) = regexprep(x(~isgood), '^([0-9\-\+\s]*)(.*)', '$2$1');
x = strtrim(x);

% After that, rely on Matlab's autoconvert

x = matlab.lang.makeValidName(x);