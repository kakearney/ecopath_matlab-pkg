function [EM, A, Eii] = eiixml2ecopathmodel(file)
%EIIXML2ECOPATHMODEL Create ecopathmodel object from EwE6 xml file
%
% EM = eiixml2ecopathmodel(file)
% [EM, A] = eiixml2ecopathmodel(file)
%
% This function creates an ecopathmodel object based on the data in an EwE6
% database file.  The .eiixml format is an alternative export option
% provided by the EwE software to save model data without reliance on
% Microsft Access drivers.  The ecopath_matlab suite does not
% support direct reading of Microsoft Access 2007 (.accdb files), so export
% through the .eiixml format is the recommended method for transfering a
% model from EwE6.6+ to the ecopathmodel format. 
%
% Input variables:
%
%   file:   name of Ecopath file.  Ecopath 6 uses Microsoft Access 2003
%           database files to store its data; by default they have .ewemdb
%           extensions
%
% Output variables:
%
%   EM:     ecopathmodel object
%
%   A:      structure of table arrays holding all data from the file in
%           its original tables.
%
%   Eii:    structure with raw xml data, as read directly from the file 
%           with Matlab's native structured data reader (readstruct)

% Copyright 2023 Kelly Kearney

%-------------------------------
% Extract data from .eiixml file
%-------------------------------

if ~exist(file, 'file')
    error('Could not find file: %s', file);
end

% Read XML

Eii = readstruct(file, 'filetype', 'xml');

A = struct;
for it = 1:length(Eii.Table)
    A = parseeiitable(Eii.Table(it), A);
end

%----------------------------
% Build ecopathmodel object
%----------------------------

% Reorder data if necessary (had to make things difficult, didn't you?).

[A.EcopathGroup, isrtg] = sortrows(A.EcopathGroup, 'Sequence');
[A.EcopathFleet, isrtf] = sortrows(A.EcopathFleet, 'Sequence');

% Group, living group, and gear/fleet number

ngroup = size(A.EcopathGroup,1);
nlive = sum(A.EcopathGroup.Type ~= 2);
ngear  = size(A.EcopathFleet,1);

% Group, fleet, and stanza names

if ~isempty(A.Stanza)
    stz = alternames(A.Stanza.StanzaName);
else
    stz = {};
end

grp = alternames(A.EcopathGroup.GroupName);
flt = alternames(A.EcopathFleet.FleetName);

old = [A.EcopathGroup.GroupName; A.EcopathFleet.FleetName; A.Stanza.StanzaName];
new = [grp; flt; stz];
issame = strcmp(old, new);
if ~all(issame)
    tmp = [old(~issame) new(~issame)]';
    str = sprintf('  %s -> %s\n', tmp{:});
    warning('Ecopathmodel:mdb2ecopathmodel:names', 'Some names changed to meet Matlab''s variable name restrictions:\n%s', str);
end

% Set up ecopathmodel object

EM = ecopathmodel(ngroup, nlive, ngear, ...
    'groups', grp, ...
    'fleets', flt, ...
    'pp', A.EcopathGroup.Type, ...
    'stanzas', stz);

% Basic info

EM.groupdata.areafrac = A.EcopathGroup.Area;
EM.groupdata.b        = A.EcopathGroup.Biomass;
EM.groupdata.pb       = A.EcopathGroup.ProdBiom;
EM.groupdata.qb       = A.EcopathGroup.ConsBiom;
EM.groupdata.ee       = A.EcopathGroup.EcoEfficiency;
EM.groupdata.ge       = A.EcopathGroup.ProdCons;
EM.groupdata.gs       = A.EcopathGroup.Unassim;
EM.groupdata.dtImp    = A.EcopathGroup.DtImports;

EM.groupdata.ee(EM.nlive + (1:(EM.ngroup-EM.nlive))) = NaN; % Correct errant 0 from file

% Imports

EM.groupdata.import = A.EcopathGroup.ImpVar;

% Diet

[tf, seqi] = ismember(A.EcopathDietComp.PreyID, A.EcopathGroup.GroupID);
[tf, seqj] = ismember(A.EcopathDietComp.PredID, A.EcopathGroup.GroupID);

dc = full(sparse(seqi, seqj, double(A.EcopathDietComp.Diet), ngroup, ngroup));
dc(:,EM.groupdata.pp>1) = 0; % Some files seem to store flow to det fractions here

EM.dc(:,:) = num2cell(dc);

% Detritus fate

detidx = A.EcopathGroup.Sequence(A.EcopathGroup.Type == 2);

df = full(sparse(seqi, seqj, double(A.EcopathDietComp.DetritusFate), ngroup, ngroup));
df = df(detidx,:)'; % Why is this backwards?... so confusing

EM.df(:,:) = num2cell(df);

% Immigration/emigration

EM.groupdata.immig    = A.EcopathGroup.Immigration;

userate = A.EcopathGroup.EmigRate ~= 0 & A.EcopathGroup.Emigration == 0;
EM.groupdata.emig(~userate)    = A.EcopathGroup.Emigration(~userate);
EM.groupdata.emigRate(userate) = A.EcopathGroup.EmigRate(userate);

% EM.groupdata.emig     = A.EcopathGroup.Emigration;
% EM.groupdata.emigRate = A.EcopathGroup.EmigRate;

% Biomass accumulation

userate = A.EcopathGroup.BiomAccRate ~= 0 & A.EcopathGroup.BiomAcc == 0;
EM.groupdata.ba(~userate)    = A.EcopathGroup.BiomAcc(~userate);
EM.groupdata.baRate(userate) = A.EcopathGroup.BiomAccRate(userate);

% EM.groupdata.ba       = A.EcopathGroup.BiomAcc;
% EM.groupdata.baRate   = A.EcopathGroup.BiomAccRate;

% Landings and discards

[tf, seqci] = ismember(A.EcopathCatch.GroupID, A.EcopathGroup.GroupID);
[tf, seqcf] = ismember(A.EcopathCatch.FleetID, A.EcopathFleet.FleetID);

landing = full(sparse(seqci, seqcf, double(A.EcopathCatch.Landing),  ngroup, ngear));
discard = full(sparse(seqci, seqcf, double(A.EcopathCatch.Discards), ngroup, ngear));

EM.landing(:,:) = num2cell(landing);
EM.discard(:,:) = num2cell(discard);

% Discard fate

[tf, seqdi] = ismember(A.EcopathDiscardFate.GroupID, A.EcopathGroup.GroupID);
[tf, seqdf] = ismember(A.EcopathDiscardFate.FleetID, A.EcopathFleet.FleetID);

discardFate = full(sparse(seqdi, seqdf, double(A.EcopathDiscardFate.DiscardFate), ngroup, ngear));
discardFate = discardFate(detidx,:)';

EM.discardFate(:,:) = num2cell(discardFate);

% Multi-stanza set details

[tf, seqs] = ismember(A.StanzaLifeStage.GroupID, A.EcopathGroup.GroupID);

EM.groupdata.stanza = zeros(ngroup,1);
EM.groupdata.stanza(seqs) = A.StanzaLifeStage.StanzaID;

EM.groupdata.ageStart = zeros(ngroup,1);
EM.groupdata.ageStart(seqs) = A.StanzaLifeStage.AgeStart;

EM.stanzadata.stanzaID = A.Stanza.StanzaID;
EM.stanzadata.BABsplit = A.Stanza.BABsplit;
EM.stanzadata.WmatWinf = A.Stanza.WmatWinf;
EM.stanzadata.RecPower = A.Stanza.RecPower;

EM.groupdata.vbK = A.EcopathGroup.vbK;
if iscell(EM.groupdata.vbK)
    % I think this only occurs when no multistanza groups present
    isemp = cellfun('isempty', EM.groupdata.vbK);
    if ~all(isemp)
        error('TODO: Found something new in vbK field... fix! (or contact Kelly to fix!)');
    end
    EM.groupdata.vbK = nan(size(EM.groupdata.vbK));
end
EM.groupdata.vbK(EM.groupdata.vbK == -1) = NaN;

% Pedigree: EwE6 allows pedigree values for B, PB, QB, DC, and catch, on a
% per group basis.  For diet, the same value is applied to all prey
% components of a group's diet.  For catch, the same value is applied to
% the part of each gear's landing and/or discard of that group.

if isfield(A, 'EcopathGroupPedigree') && ~isempty(A.EcopathGroupPedigree)

    tbl = {...
        'BiomassAreaInput'  'b'   
        'Biomass'           'b'    % ?? Found in Albatross Bay, upgraded from older EwE model
        'PBInput'           'pb'    
        'QBInput'           'qb'    
        'DietComp'          'dc'
        'TCatchInput'       'catch'};

    [~,loc] = ismember(A.EcopathGroupPedigree.VarName, tbl);
    vped = tbl(loc,2);

    [~, grpidx] = ismember(A.EcopathGroupPedigree.GroupID, A.EcopathGroup.GroupID);
    [~, lidx] = ismember(A.EcopathGroupPedigree.LevelID, A.Pedigree.LevelID);
    lped = A.Pedigree.Confidence(lidx)./100;    

    [isgd, cidxgd] = ismember(vped, EM.groupdata.Properties.VariableNames);

    Ped.property = repmat({'groupdata'}, sum(isgd), 1);
    Ped.row      = grpidx(isgd);
    Ped.column   = cidxgd(isgd);
    Ped.pedigree = lped(isgd);

    isdc = strcmp(vped, 'dc');

    dcped = kron(lped(isdc), ones(ngroup,1));
    dcrow = repmat((1:ngroup)', sum(isdc), 1);
    dccol = kron(grpidx(isdc), ones(ngroup,1));

    idx = sub2ind(size(dc), dcrow, dccol);  
    mask = dc(idx) > 0;

    Ped.property = [Ped.property; repmat({'dc'}, sum(mask), 1)];
    Ped.row      = [Ped.row;      dcrow(mask)];
    Ped.column   = [Ped.column;   dccol(mask)];
    Ped.pedigree = [Ped.pedigree; dcped(mask)];

    isfsh = strcmp(vped, 'catch');
    fshcatch = landing + discard;

    fshped = kron(lped(isfsh), ones(ngear,1));
    fshrow = kron(grpidx(isfsh), ones(ngear,1));
    fshcol = repmat((1:ngear)', sum(isfsh), 1);

    idx = sub2ind(size(fshcatch), fshrow, fshcol);
    mask = landing(idx) > 0;

    Ped.property = [Ped.property; repmat({'landing'}, sum(mask), 1)];
    Ped.row      = [Ped.row;      fshrow(mask)];
    Ped.column   = [Ped.column;   fshcol(mask)];
    Ped.pedigree = [Ped.pedigree; fshped(mask)];

    mask = discard(idx) > 0;

    Ped.property = [Ped.property; repmat({'discard'}, sum(mask), 1)];
    Ped.row      = [Ped.row;      fshrow(mask)];
    Ped.column   = [Ped.column;   fshcol(mask)];
    Ped.pedigree = [Ped.pedigree; fshped(mask)];

    EM.pedigree = struct2table(Ped);
end


%----------------------------
% Variable name adjustments
%----------------------------

function x = alternames(x)

isgood = cellfun(@isvarname, x);

% Often named using ages, so move numeric bits to end

x(~isgood) = regexprep(x(~isgood), '^([0-9\-\+\s]*)(.*)', '$2$1');
x = strtrim(x);

% After that, rely on Matlab's autoconvert

x = matlab.lang.makeValidName(x);

%----------------------------
% Parse xml tables
%----------------------------

function A = parseeiitable(T, A)

cols = reshape(strsplit(T.ColumnsAttribute, {',', ':'}), 2, []);
ncol = size(cols,2);
if ismissing(T.Row)
    nrow = 0;
else
    nrow = length(T.Row);
end

if nrow > 0
    fmtlookup = {...
        'System.Single'     '%f32'
        'System.Int64'      '%d64'
        'System.Int32'      '%d32' 
        'System.Byte'       '%bu8'
        'System.Boolean'    '%s' % <- needs conversion after
        'System.String'     '%s'
        };
    
    [tf,loc] = ismember(cols(2,:), fmtlookup(:,1));
    if ~all(tf)
        error('TODO: need to add new format to lookup table!');
    end
    fmtstr = sprintf('%s', fmtlookup{loc,2});
    
    if nrow == 1 && loc(end)==6
        str = T.Row + ","; % allows for empty last string
    else
        str = strjoin(T.Row, ",");
    end
    str = strrep(str, 'Infinity', 'Inf');

    data = textscan(str, fmtstr, 'delimiter', ',');
    if any(loc == 5) % logical
        data(loc==5) = cellfun(@(x) strcmp(x,'True'), data(loc==5), 'uni', 0);
    end
    
    A.(T.NameAttribute) = table(data{:}, 'variablenames', cols(1,:));
else
    A.(T.NameAttribute) = cell2table(cell(0,ncol), 'variablenames', cols(1,:));
end



