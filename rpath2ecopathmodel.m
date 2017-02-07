function [EM,B] = rpath2ecopathmodel(varargin)
%RPATH2ECOPATHMODEL Create ecopathmodel object from Rpath data files
%
% [A, B] = rpath2ecopathmodel(modfile, dietfile)
% [A, B] = rpath2ecopathmodel(modfile, dietfile, p1, v1, ...)
%
% Rpath (the R-based version of Ecopath/Ecosim, available from
% https://github.com/slucey/RpathDev) allows model data be imported from
% and exported to .csv files.  This function reads data from a set of these
% files into an ecopathmodel object.
%
% A model *must* include at least two files: the main model table (see
% modfile input) and the diet table (see dietfile input).  All other files
% are optional (though the stanza and stanza group tables should always be
% present or absent as a pair).  If those files are empty, this function
% assumes that they are not relevant to the model being imported. Note that
% the write.rpath.params R function in the Rpath package exports all
% tables, regardless of whether they are empty or not; this function can
% handle those empty tables if you choose to list them as input parameters,
% but they aren't necessary. 
% 
% Input variables:
%
%   modfile:            Name of base model data table file (full path, with
%                       .csv extension)  
%                       columns: 
%
%   dietfile:           Name of diet data table file (full path, with .csv
%                       extension)   
%
% Optional input variables (passed as parameter/value pairs)
%
%   pedfile:            Name of pedigree data table file (full path, with
%                       .csv extension).  If empty, no pedigree data will
%                       be read. Please note that the pedigree values used
%                       by Rpath right now (by default, a table of all
%                       ones) are just placeholders, and are unlikely to be
%                       appropriate to an ecopathmodel.  This may change in
%                       the future as Rpath development continues, but for
%                       now I recommend not using that table at all when
%                       importing Rpath data.
%                       [''] 
%
%   stanzafile:         Name of stanza table file (full path, with .csv
%                       extension).  Must be accompanied by a
%                       stanzagroupfile input.  If empty, assumes that no
%                       groups refer to multi-stanza sets. ['']
%
%   stanzagroupfile:    Name of stanza groups table file (full path, with
%                       .csv extension).  Must be accompanied by a
%                       stanzafile input.  If empty, assumes that no groups
%                       refer to multi-stanza sets. ['']
%
%   juvsfile:           Name of old-style adult/juveniles table file (full
%                       path, with .csv extension).  This type of file was
%                       used by older versions of Rpath, prior to the
%                       introduction of stanza/stanzagroup tables, so this
%                       code can also use it in place of the
%                       stanza/stanzagroup pair of files. [''] 
%
% Output variables:
%
%   EM:         ecopathmodel object
%
%   B:          structure with full input tables from files

% Copyright 2016-2017 Kelly Kearney

%--------------------
% Parse input
%--------------------

p = inputParser;
p.addRequired( 'modfile',             @(x) validateattributes(x, {'char'}, {}));
p.addRequired( 'dietfile',            @(x) validateattributes(x, {'char'}, {}));
p.addParameter('pedfile',         '', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('stanzagroupfile', '', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('stanzafile',      '', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('juvsfile',        '', @(x) validateattributes(x, {'char'}, {}));
p.parse(varargin{:});

Opt = p.Results;

% Check that files exist

if ~exist(Opt.modfile, 'file')
    error('Could not find specified model file (%s); check name', Opt.modfile);
end

if ~exist(Opt.dietfile, 'file')
    error('Could not find specified diet file (%s); check name', Opt.dietfile);
end

if ~isempty(Opt.pedfile) && ~exist(Opt.pedfile, 'file')
    error('Could not find specified pedigree file (%s); check name', Opt.pedfile);
end

if ~isempty(Opt.stanzagroupfile) && ~exist(Opt.stanzagroupfile, 'file')
    error('Could not find specified stanza group file (%s); check name', Opt.stanzagroupfile);
end

if ~isempty(Opt.stanzafile) && ~exist(Opt.stanzafile, 'file')
    error('Could not find specified stanza file (%s); check name', Opt.stanzafile);
end

% Check that stanza files as present/absent together

if (~isempty(Opt.stanzafile) &&  isempty(Opt.stanzafile)) || ...
   ( isempty(Opt.stanzafile) && ~isempty(Opt.stanzafile))    
    error('Stanza and stanza group files must be provided as a pair');
end

if ~isempty(Opt.juvsfile)
    if ~isempty(Opt.stanzafile) || ~isempty(Opt.stanzagroupfile)
        error('Multi-stanza groups can be defined by a juvs file or stanza+stanzagroup files, not both');
    end
    Opt.old = true;
else
    Opt.old = false;
end

hasped = ~isempty(Opt.pedfile);
hasstz = ~isempty(Opt.stanzafile) || ~isempty(Opt.juvsfile);

% p = inputParser;
% p.addParameter('basestr', '_model',    @(x) validateattributes(x, {'char'}, {}));
% p.addParameter('dietstr', '_diet',     @(x) validateattributes(x, {'char'}, {}));
% p.addParameter('juvsstr', '_juvs',     @(x) validateattributes(x, {'char'}, {}));
% p.addParameter('pedstr',  '_pedigree', @(x) validateattributes(x, {'char'}, {}));
% p.addParameter('stanzastr', '_stanzas', @(x) validateattributes(x, {'char'}, {}));
% p.addParameter('stgrpstr', '_stanza_groups', @(x) validateattributes(x, {'char'}, {}));
% p.addParameter('old', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
% p.parse(varargin{:});
% 
% Opt = p.Results;

%--------------------
% Read tables
%--------------------

% Read file data into tables

Base = readtable(Opt.modfile, 'TreatAsEmpty', 'NA');
Diet = readtable(Opt.dietfile, 'ReadRowNames', true, 'TreatAsEmpty', 'NA');

if hasped
    Ped  = readtable(Opt.pedfile, 'TreatAsEmpty', 'NA');
end

if hasstz
    if Opt.old
        Juvs = readtable(Opt.juvsfile);

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
        Stan = readtable(Opt.stanzafile, 'TreatAsEmpty', 'NA');
        Sgrp = readtable(Opt.stanzagroupfile, 'TreatAsEmpty', 'NA');

        isn = isnan(Stan.StGroupNum);
        Stan = Stan(~isn,:);

        isn = isnan(Sgrp.StGroupNum);
        Sgrp = Sgrp(~isn,:);
    end
end

B.Base = Base;
B.Diet = Diet;
if hasped
    B.Ped = Ped;
else
    B.Ped = [];
end
if hasstz
    B.Stan = Stan;
    B.Sgrp = Sgrp;
    if Opt.old
        B.Juvs = Juvs;
    else
        B.Juvs = [];
    end
else
    B.Stan = [];
    B.Sgrp = [];
    B.Juvs = [];
end

% sflag = ~isempty(Stan);

% Parse column names

bcol = parsecolname(Base);
dcol = parsecolname(Diet);
if hasped
    pcol = parsecolname(Ped);
end

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
if hasstz
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

if hasstz
    EM.groupdata.stanza(Stan.GroupNum) = Stan.StGroupNum;
    EM.groupdata.ageStart(Stan.GroupNum) = Stan.First;

    [tf, loc] = ismember(EM.groupdata.stanza, Sgrp.StGroupNum);

    EM.groupdata.vbK(tf) = Sgrp.VBGF_Ksp(loc(tf));

    EM.stanzadata.stanzaID = Sgrp.StGroupNum;
end

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

if hasped
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
end

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