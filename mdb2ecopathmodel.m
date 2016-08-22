function [EM, A] = mdb2ecopathmodel(file)
%MDB2ECOPATHMODEL Create ecopathmodel object from EwE6 data file
%
% EM = mdb2ecopathmodel(file)
% [EM, A] = mdb2ecopathmodel(file)
%
% This function creates an ecopathmodel object based on the data in an EwE6
% database file.  It is only designed to work with EwE version 6 files (has
% been tested for EwE6.4 and EwE6.5); older files must be converted via the
% conversion utility that comes with EwE6.  The newer files are usually
% saved with a .ewemdb (or .EwEmdb) file extension, as opposed to older
% ones that simply carried the .mdb extension.
%
% This function relies on the mdbtools utilities to read the MS Access
% database files (https://github.com/brianb/mdbtools).  Mac users can get
% it via either MacPorts or Homebrew if they don't want to compile from
% source. Please make sure this utility is properly compiled prior to
% calling mdb2ecopathmodel.m.  Windows users: you will need a C compiler (I
% believe you can get one for free via Visual Studio Express).
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
%   A:      structure of dataset arrays holding all data from the file in
%           its original tables.

% Copyright 2016 Kelly Kearney

%----------------------------
% Extract data from .mdb file
%----------------------------

if ~exist(file, 'file')
    error('Could not find file: %s', file);
end

% Translate file names such that they can be passed via command line

file = regexprep(file, '([\s,<>|:\(\)&;\?\*])', '\\$1');

% Read table names

cmd = sprintf('mdb-tables -d, %s', file);
[s,r] = system(cmd);
if s
    if s == 127
        msg = 'Error calling mdb-tables.  Please make sure you have installed the mdbtools utility and that it is accessible to Matlab (Note that your MATLAB system path, seen with getenv(''PATH''), may differ from your Matlab search path; the former is what matters here)';
        msg = wraptext(msg, 70);
        error('%s', msg);
    else
        error('Error reading table names: %s', r);
    end
end
    
tables = regexp(r, ',', 'split');
tables = regexprep(tables, '\n', '');
isemp = cellfun('isempty', tables);
tables = tables(~isemp);  

% Crazy regular expression for Excel-style comma-delimited stuff... find
% commas that aren't imbedded within quotes.  No longer needed with
% readtable and readtext... but keeping it for reference because I'll never
% be able to reinvent this ridiculousness. :-)

pattern = ',(?=(?:[^\"]*\"[^\"]*\")*(?![^\"]*\"))';

% Read all Ecopath table data

err = cell(0,2);
for it = 1:length(tables)
    
    mdbtmp = [tempname '.txt'];
    
    if regexpfound(tables{it}, '\s')
        tbl = regexprep(tables{it}, '([\s,<>|:\(\)&;\?\*])', '\\$1');
        cmd = sprintf('mdb-export %s %s > %s', file, tbl, mdbtmp);
        tables{it} = regexprep(tables{it}, '\s', '_');
    else
        cmd = sprintf('mdb-export %s %s > %s', file, tables{it}, mdbtmp);
    end
    [s,r] = system(cmd);
    if s
        err = [err; {tables{it} r}];
        continue
    end
    
    try
        
        A.(tables{it}) = readtable(mdbtmp, 'Delimiter', ',');
        vname = A.(tables{it}).Properties.VariableNames;
        
        for iv = 1:length(vname)
            if isnumeric(A.(tables{it}).(vname{iv}))
                A.(tables{it}).(vname{iv})(A.(tables{it}).(vname{iv}) == -9999) = NaN;
            else
                A.(tables{it}).(vname{iv}) = regexprep(A.(tables{it}).(vname{iv}), {'(^")|("$)', '""'}, {'', '"'}); % Strip out beginning/end quotes and escaped quotes
            end
        end
        
    catch
        
        try
            % readtable will get tripped up if a field includes a newline
            % character enclosed within quotes, so we need to manually
            % parse that.

            [data, s] = readtext(mdbtmp, ',', '', '"');
            data(s.stringMask) = regexprep(data(s.stringMask), {'(^")|("$)', '""'}, {'', '"'});
            
            tmp = cell2mat(data(s.numberMask));
            tmp(tmp == -9999) = NaN;
            data(s.numberMask) = num2cell(tmp);
            
            A.(tables{it})  = cell2table(data(2:end,:), 'VariableNames', data(1,:));
            
        catch
            
            % if still no good, just return the text (and warn if the user
            % is getting this info)
            
            if nargout > 1
                warning('Ecopathmodel:mdb2ecopathmodel:parse', 'Could not parse table: %s', tables{it});
            end
            A.(tables{it}) = r;
            
        end
    end
end

if ~isfield(A, 'EcopathGroup')
    error('Expected tables not found in file; please check that this is an EwE6-formatted file');
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

dc = full(sparse(seqi, seqj, A.EcopathDietComp.Diet, ngroup, ngroup));
dc(:,EM.groupdata.pp>1) = 0; % Some files seem to store flow to det fractions here

EM.dc(:,:) = num2cell(dc);

% Detritus fate

detidx = A.EcopathGroup.Sequence(A.EcopathGroup.Type == 2);

df = full(sparse(seqi, seqj, A.EcopathDietComp.DetritusFate, ngroup, ngroup));
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

landing = full(sparse(seqci, seqcf, A.EcopathCatch.Landing,  ngroup, ngear));
discard = full(sparse(seqci, seqcf, A.EcopathCatch.Discards, ngroup, ngear));

EM.landing(:,:) = num2cell(landing);
EM.discard(:,:) = num2cell(discard);

% Discard fate

[tf, seqdi] = ismember(A.EcopathDiscardFate.GroupID, A.EcopathGroup.GroupID);
[tf, seqdf] = ismember(A.EcopathDiscardFate.FleetID, A.EcopathFleet.FleetID);

discardFate = full(sparse(seqdi, seqdf, A.EcopathDiscardFate.DiscardFate, ngroup, ngear));
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

if ~isempty(A.EcopathGroupPedigree)

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


