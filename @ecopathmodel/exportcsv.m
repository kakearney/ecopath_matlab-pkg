function exportcsv(EM, filename)
%EXPORTCSV Export ecopathmodel object data to csv files
%
% exportcsv(EM, filename)
%
% This function exports the data in an ecopathmodel object to a set of
% comma-delimited text files.  The files are formatted to match the
% write.rpath.params method from the Rpath R package, allowing models to be
% used interchangeably.  
%
% Input variables:
%
%   EM:     ecopathmodel object
%   
%   fname:  base name for output files. Four files will be produced:
%           fname_model.csv:            Includes most group-related
%                                       parameters  
%           fname_diet.csv:             Includes diet fraction matrix
%           fname_stanzas.csv:          Includes growth curve parameters
%                                       for stanza groups 
%           fname_stanza_groups.csv:    Includes growth curve parameters
%                                       for stanza sets 

% Copyright 2016 Kelly Kearney


% Check input

narginchk(2,2);

validateattributes(filename, {'char'}, {});

[pth,fl,ex] = fileparts(filename);
if ~isempty(pth) && ~exist(pth, 'dir')
    error('The folder %s does not exist', pth);
end
if ~isempty(ex) && strcmp(ex, '.csv')
    filename = fullfile(pth, fl);
end

% Basic parameters

colnames = {'Group', 'Type', 'Biomass', 'PB', 'QB', 'EE', 'ProdCons', ...
    'BioAcc', 'Unassim', 'DetInput'};
nbasic = length(colnames)-1;

isdet = EM.groupdata.pp == 2;
ndet = EM.ngroup - EM.nlive;

discname = cellfun(@(x) [x '.disc'], EM.fleet, 'uni', 0);

colnames = [colnames EM.name(isdet)' EM.fleet' discname'];
rownames = [EM.name; EM.fleet];

modeldata = nan(EM.ngroup+EM.ngear, nbasic+ndet+EM.ngear*2);

modeldata(1:EM.ngroup,1) = EM.groupdata.pp;
modeldata((1:EM.ngear)+EM.ngroup,1) = 3;

modeldata(1:EM.ngroup,2) = EM.groupdata.b;
modeldata(1:EM.ngroup,3) = EM.groupdata.pb;
modeldata(1:EM.ngroup,4) = EM.groupdata.qb;
modeldata(1:EM.ngroup,5) = EM.groupdata.ee;
modeldata(1:EM.ngroup,6) = EM.groupdata.ge;
modeldata(1:EM.ngroup,7) = EM.groupdata.ba;
modeldata(1:EM.ngroup,8) = EM.groupdata.gs;

modeldata((EM.nlive+1):EM.ngroup,9) = EM.groupdata.dtImp((EM.nlive+1):EM.ngroup);

modeldata(:,nbasic+(1:ndet)) = [table2array(EM.df); table2array(EM.discardFate)];

modeldata(1:EM.ngroup,(nbasic+ndet+1):end) = [table2array(EM.landing) table2array(EM.discard)];

isn = isnan(modeldata);
modeldata = num2cell(modeldata);
[modeldata{isn}] = deal('NA');

modeldata = [colnames; [rownames modeldata]];
isnum = cellfun(@isnumeric, modeldata);
modeldata(isnum) = cellfun(@num2str, modeldata(isnum), 'uni', 0);

% Diet

dc = table2array(EM.dc);
is0 = dc == 0;
dietdata = num2cell(dc);
[dietdata{is0}] = deal('NA');
dietdata(~is0) = cellfun(@num2str, dietdata(~is0), 'uni', 0);

dietdata = [['Group' EM.name']; [EM.name dietdata]];

% Stanzas

idx = EM.stanzaindices;
stanzadata = {'StGroupNum', 'Stanza', 'GroupNum', 'Group', 'First', 'Last', 'Z', 'Leading'};
for is = 1:length(idx)
    n = length(idx{is});
    stgrpnum = ones(n,1)*is;
    stz = (1:n)';
    grpnum = idx{is};
    grp = EM.name(idx{is});
    first = EM.groupdata.ageStart(idx{is});
    last = [first(2:end)-1; 400];
    z = EM.groupdata.pb(idx{is});
    leading = cell(n,1);
    leading{end} = 'TRUE';
    [leading{1:end-1}] = deal('FALSE');
    
    stanzadata = [stanzadata; ...
        [num2cell([stgrpnum stz grpnum]) grp num2cell([first last z]) leading]];
    
end

isnum = cellfun(@isnumeric, stanzadata);
stanzadata(isnum) = cellfun(@num2str, stanzadata(isnum), 'uni', 0);

% Stanza groups

stgrpcol = {'StGrpNum', 'StanzaGroup', 'nstanzas', 'VGBF_Ksp', 'VGBF_d', 'Wmat', 'RecPower'};

stgrpnum = EM.stanzadata.stanzaID;
stgrp = EM.stanza;
nst = cellfun(@length, idx);
vbK = cellfun(@(x) EM.groupdata.vbK(x(1)), idx);
vbd = ones(size(vbK)) * (2/3); % ??
wmat = EM.stanzadata.WmatWinf;
recpower = EM.stanzadata.RecPower;

stgrpdata = [stgrpcol; [num2cell(stgrpnum) stgrp num2cell([nst vbK vbd wmat recpower])]];
isnum = cellfun(@isnumeric, stgrpdata);
stgrpdata(isnum) = cellfun(@num2str, stgrpdata(isnum), 'uni', 0);

% Write to files

modelfile  = [filename '_model.csv'];
dietfile   = [filename '_diet.csv'];
stanzafile = [filename '_stanzas.csv'];
stgrpfile  = [filename '_stanza_groups.csv'];

fmt = repmat('%s,', 1, nbasic+1+ndet+EM.ngear*2);
modeldata = modeldata';
fid = fopen(modelfile, 'w+');
fprintf(fid, [fmt(1:end-1) '\n'], modeldata{:});
fclose(fid);

fmt = repmat('%s,', 1, EM.ngroup+1);
dietdata = dietdata';
fid = fopen(dietfile, 'w+');
fprintf(fid, [fmt(1:end-1) '\n'], dietdata{:});
fclose(fid);

fmt = repmat('%s,', 1, size(stanzadata,2));
stanzadata = stanzadata';
fid = fopen(stanzafile, 'w+');
fprintf(fid, [fmt(1:end-1) '\n'], stanzadata{:});
fclose(fid);

fmt = repmat('%s,', 1, size(stgrpdata,2));
stgrpdata = stgrpdata';
fid = fopen(stgrpfile, 'w+');
fprintf(fid, [fmt(1:end-1) '\n'], stgrpdata{:});
fclose(fid);


