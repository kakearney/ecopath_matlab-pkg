function ecopathmodel2rpath(EM, outbase)
%ECOPATHMODEL2RPATH Print ecopathmodel data to comma-delimited files
%
% ecopathmodel2rpath(EM, outbase)
%
% This function saves ecopathmodel data to a set of comma-delimited files,
% formatted as required by the Rpath R package.  These files can be read
% into Matlab via rpath2ecopathmodel, or into R via read.rpath.param.
%
% Input variables:
%
%   EM:         ecopathmodel object
%
%   outbase:    base name for output files (can include full path, but
%               folder must exist prior to calling this function).  Files
%               will be saved with this name, plus the extensions
%               '_model.csv' and '_diet.csv' (and '_stanzas.csv', and
%               '_stanza_groups.csv', if multi-stanza sets are present).
%
%               Note: at this time, the Rpath pedigree table is only a
%               placeholder, and the values in it are never used by the
%               rpath or rsim R functions.  Therefore, I do not export
%               ecopathmodel pedigree values to file.  This may change in
%               the future as Rpath development continues.

% Copyright 2016 Kelly Kearney

%----------------------
% Setup
%----------------------

isdet = EM.groupdata.pp == 2;
hasstz = ~isempty(EM.stanza);

addquotes = @(x) cellfun(@(a) sprintf('"%s"',a),x, 'uni', 0); 

% Check if path exists

pth = fileparts(outbase);
if ~isempty(pth) && ~exist(pth, 'dir')
    error('Output folder not found');
end

%----------------------
% Data tables
%----------------------

% Model data table

modeldata = [[EM.groupdata.pp; ones(EM.ngear,1)*3] ...
             [EM.groupdata.b; nan(EM.ngear,1)] ...
             [EM.groupdata.pb; nan(EM.ngear,1)] ...
             [EM.groupdata.qb; nan(EM.ngear,1)] ...
             [EM.groupdata.ee; nan(EM.ngear,1)] ...
             [EM.groupdata.ge; nan(EM.ngear,1)] ...
             [EM.groupdata.ba; nan(EM.ngear,1)] ...
             [EM.groupdata.gs; nan(EM.ngear,1)] ...
             [EM.groupdata.dtImp; nan(EM.ngear,1)] ...
             [table2array(EM.df); table2array(EM.discardFate)] ...
             [table2array(EM.landing); nan(EM.ngear)] ...
             [table2array(EM.discard); nan(EM.ngear)]];
         
colname = addquotes(...
    [{'Group', 'Type', 'Biomass', 'PB', 'QB', 'EE', 'ProdCons', 'BioAcc', 'Unassim', 'DetInput'}, ...
    EM.name(isdet)', ...
    EM.fleet', ...
    cellfun(@(x) [x '.disc'], EM.fleet', 'uni', 0)]);

rowname = addquotes([EM.name; EM.fleet]);

modeldata = [colname; [rowname num2cell(modeldata)]];

% Diet data table  

colname = addquotes(['Group'; EM.name(~isdet)]');
rowname = addquotes([EM.name; 'Import']);

dietdata = [table2array(EM.dc(:,~isdet)); EM.groupdata.import(~isdet)'];
dietdata = [colname; [rowname num2cell(dietdata)]];

if hasstz

    % Stanza group data table

    colname = addquotes({'StGroupNum', 'StanzaGroup', 'nstanzas', 'VGBF_Ksp', 'VGBF_d', 'Wmat', 'RecPower'});

    sid = EM.stanzaindices;

    nstz = cellfun(@length, sid);
    s1 = cellfun(@(x) x(1), sid);

    stgrpdata = [colname; [num2cell(EM.stanzadata.stanzaID), EM.stanza, ...
        num2cell([nstz EM.groupdata.vbK(s1) ones(length(nstz),1)*(2/3), ...
        EM.stanzadata.WmatWinf, EM.stanzadata.RecPower])]];

    % Stanzas data table

    colname = addquotes({'StGroupNum', 'StanzaNum', 'GroupNum', 'Group', 'First', 'Last', 'Z', 'Leading'});

    sgnum = cell2mat(arrayfun(@(n,x) ones(n,1)*x, nstz, (1:length(sid))', 'uni', 0));
    snum = cell2mat(arrayfun(@(x) (1:x)', nstz, 'uni', 0));

    sid = cell2mat(sid);

    age1 = EM.groupdata.ageStart(sid);
    maxage = 1200;
    age2 = [age1(2:end)-1; maxage]; % Arbitrary choice of 100 years for upper bound
    islead = age2 == -1;
    age2(islead) = maxage;

    stanzadata = [colname; ...
        [num2cell([sgnum snum sid]), ...
         EM.name(sid), ...
         num2cell([age1 age2 EM.groupdata.pb(sid) islead]), ...
        ]];
end

%----------------------
% Save to files
%----------------------

printcsv( modeldata, [outbase '_model.csv']);
printcsv(  dietdata, [outbase '_diet.csv']);
if hasstz
    printcsv(stanzadata, [outbase '_stanzas.csv']);
    printcsv( stgrpdata, [outbase '_stanza_groups.csv']);
end

%----------------------
% Subfunction: format
% and save
%----------------------

function x = printcsv(x, file)

% Convert all cell contents to strings

islog = cellfun(@islogical, x);
isna = cellfun(@(x) isnumeric(x) && isnan(x), x);
isnum = cellfun(@isnumeric, x);

tf = {'FALSE', 'TRUE'};

x(islog) = cellfun(@(a) tf(a+1), x(islog), 'uni', 0);
[x{isna}] = deal('NA');
x(isnum) = cellfun(@(a) num2str(a), x(isnum), 'uni', 0);

% Print as comma-delimited text

ncol = size(x,2);
fmt = repmat('%s,', 1, ncol);
fmt = [fmt(1:end-1) '\n'];

x = x';

fid = fopen(file, 'wt');
fprintf(fid, fmt, x{:});
fclose(fid);
