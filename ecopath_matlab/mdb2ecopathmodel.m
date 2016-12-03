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
% In order to extract data from an MS Access database, this function needs
% to use an ODBC database driver.  It does so in one of two ways, depending
% on your operating system:
%
% - Linux/MacOS: mdbtools utilities
% - Windows: pyodbc python module, using locally-installed ODBC driver
%
% See below for installation instructions for either of these options.
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
% External software installation instructions:
%
% * Linux:
% 
%   The mdbtools source code can be acquired from GitHub:
%   https://github.com/brianb/mdbtools.  Instructions to compile from
%   source are included in the download.
%
% * MacOS:
%
%   Ports of mdbtools are available through either Homebrew or MacPorts. 
%   Install either Homebrew or MacPorts if you don't already have them on
%   your system, then run:
%     brew install mdbtools       <-- from Terminal, Homebrew
%     port install mdbtools       <-- from Terminal, MacPorts
%
% * Windows:
%
%   On Windows, first check your computer to make sure you have the MS
%   Access driver installed (the driver list is buried in the Control
%   Panel... search ODBC Data Sources as a quick way to find it).  You need
%   the 'Microsoft Access Driver (*.mdb)'.  Next, you need python 3.x.
%   There are many options for downloading and installing python if you
%   don't already have it; I personally like the Anaconda Python
%   distribution (https://www.continuum.io/downloads) for scientific use.
%   Finally, download and install the pyodbc module; easy installation can
%   be done via pip:     
%     pip install pyodbc   <-- from Command Prompt 
%   or via conda, if you use Anaconda Python
%     conda install pyodbc <-- from Command Prompt 
%   My mdbexport.py script also requires the csv and os modules; these are
%   installed automatically with most commom python distributions.


% Copyright 2016 Kelly Kearney

%----------------------------
% Extract data from .mdb file
%----------------------------

if ~exist(file, 'file')
    error('Could not find file: %s', file);
end

% Check that mdb-querying tools are ready to go

checkdeps;

% Dump tables to comma-delimited files

csvfolder = tempname;
if ~exist(csvfolder, 'dir')
    mkdir(csvfolder);
end

if ispc
    
    py.mdbexport.mdbexport(csvfolder, file);
    
    tables = dir(fullfile(csvfolder, '*.csv'));
    tables = strrep({tables.name}, '.csv', '');
else
    % Translate file name such that it can be passed via command line

    file = regexprep(file, '([\s,<>|:\(\)&;\?\*])', '\\$1'); 

    % Read table names

    cmd = sprintf('mdb-tables -1 %s', file); % one per line to account for commas, spaces in names
    [s,r] = system(cmd);
    if s
        error('Error reading table names: %s', r);
    end

    tables = regexp(r, '\n', 'split');
    isemp = cellfun('isempty', tables);
    tables = tables(~isemp); 
    
    % Extract to csv
       
    for it = 1:length(tables)
        fname = fullfile(csvfolder, tables{it});
        if regexpfound(tables{it}, '\s')
            tbl = regexprep(tables{it}, '([\s,<>|:\(\)&;\?\*])', '\\$1');
            tables{it} = regexprep(tables{it}, '\s', '_');
            fname = regexprep(fname, '\s', '_');
            cmd = sprintf('mdb-export %s %s > %s.csv', file, tbl, fname);
        else
            cmd = sprintf('mdb-export %s %s > %s.csv', file, tables{it}, fname);
        end
        [s,r] = system(cmd);
        if s
            err = [err; {tables{it} r}];
        end
    end
end

% Crazy regular expression for Excel-style comma-delimited stuff... find
% commas that aren't imbedded within quotes.  No longer needed with
% readtable and readtext... but keeping it for reference because I'll never
% be able to reinvent this ridiculousness. :-)

pattern = ',(?=(?:[^\"]*\"[^\"]*\")*(?![^\"]*\"))';

% Check that this is an EwE6, not EwE5 (or any other random .mdb) file

tbls = {'EcopathGroup', 'EcopathFleet', 'Stanza', 'EcopathDietComp', ...
        'EcopathCatch', 'EcopathDiscardFate', 'StanzaLifeStage', ...
        'EcopathGroupPedigree', 'Pedigree'};
hastbl = ismember(tbls, tables);
if ~all(hastbl)
    str = sprintf('%s, ', tbls{~hastbl});
    error('Expected tables (%s) not found in file; please check that this is an EwE6-formatted file', str(1:end-1));
end

% Read all Ecopath table data

err = cell(0,2);
w = warning('off', 'MATLAB:table:ModifiedVarnames');
for it = 1:length(tables)
    mdbtmp = fullfile(csvfolder, [tables{it} '.csv']);
    try
   
        A.(tables{it}) = readtable(mdbtmp, 'Delimiter', ',');
        vname = A.(tables{it}).Properties.VariableNames;
        
        % Replace -9999 placeholders with NaNs, and strip unnecessary
        % quotes
        
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
            A.(tables{it}) = fileread(mdbtmp);
            
        end
    end
end
warning(w);

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

%----------------------------
% Dependency check
%----------------------------

function checkdeps
%CHECKDEPS Check for necessary external software
% On Linux or Mac: need mdbtools, compiled and on Matlab system path
% On Windows: need python (3+), plus pyodbc and csv modules, plus mdbexport
%             installed locally

if ispc % Windows
    [v,~,loaded] = pyversion;
    if isempty(v)
        error('Could not access python; please make sure python 3 is installed on your computer');
    end
    if str2double(v) < 3
        if loaded
            msg = 'The mdbexport.py module was written using python 3 syntax; you currently have python 2 loaded. To change versions, install python 3 if necessary, then restart Matlab and then call pyversion';
            msg = wraptext(msg, 70);
            error('%s', msg);
        else
            try
                pyversion 3.5
            catch
                try
                    pyversion 3.4
                catch
                    msg = 'The mdbexport.py module was written using python 3 syntax; I could only find version 2 on this computer.  Please install python 3 to use this function';
                    msg = wraptext(msg, 70);
                    error('%s', msg);
                end
            end
        end
    end
    modnotfound = @(ME) ~isempty(strfind(ME.message, 'No module named'));
    try
        py.importlib.import_module('pyodbc');
    catch ME
        if modnotfound(ME)
            error('Could not import pyodbc python module; please make sure it is installed');
        else
            rethrow(ME);
        end
    end
    try
        py.importlib.import_module('csv');
    catch ME
        if modnotfound(ME)
            error('Could not import csv python module; please make sure it is installed');
        else
            rethrow(ME);
        end
    end
    try
        py.importlib.import_module('mdbexport');
    catch ME
        if modnotfound(ME)
            pth = fileparts(which('mdb2ecopathmodel'));
            mdbpath = fullfile(pth, 'mdbexport', 'mdbexport');
            P = py.sys.path;
            if count(P,mdbpath) == 0
                insert(P,int32(0),mdbpath);
            end
            try
                py.importlib.import_module('mdbexport');
            catch
                error('Could not import mdbexport python module; please check that you either installed it locally or kept it in its default location in the ecopath_matlab package');
            end
        else
            rethrow(ME);
        end
    end
else    % Linux, MacOS
    [s,r] = system('mdb-ver -M');
    if s
        if s == 127
            msg = 'Could not call mdbtools function.  Please make sure you have installed the mdbtools utilities and that they are accessible to Matlab (Note that your MATLAB system path, seen with getenv(''PATH''), may differ from your Matlab search path; the former is what matters here)';
            msg = wraptext(msg, 70);
            error('%s', msg);
        else
            error('Error calling mdb-ver: %s', r);
        end
    end
end



