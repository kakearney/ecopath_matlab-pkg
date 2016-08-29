function EM = ewef2ecopathmodel(basename, varargin)
%EWEF2ECOPATHMODEL Create ecopathmodel object from EwE-F input files
%
% [A, B] = rpath2ecopathmodel(basename, p1, v1, ...)
%
% EwE-F (the Fortran-based version of Ecopath/Ecosim, available from
% https://bitbucket.org/ewe-f/) uses a set of tab-delimited text files for
% input.  This function reads data from a set of these files into an
% ecopathmodel object.  See the EwE-F User's Guide for more detail on the
% file format.
%
% Input variables:
%
%   basename:   Base file name for the input files
%
% Optional input variables (passed as paramete/value pairs)
%
%   mainstr:    string appended to basename for main parameter data file,
%               not including.txt extension ['_Scenario'] 
%
%   dietstr:    string appended to basename for diet composition parameter
%               data file, not including.txt extension ['_DC']
%
%   detfstr:    string appended to basename for detritus fate parameter
%               data file, not including.txt extension ['_DetFate']
%
%   growstr:    string appended to basename for multi-stanza group growth
%               parameter data file, not including.txt extension
%               ['_GrowthParameters']  

% Copyright 2016 Kelly Kearney

%--------------------
% Parse input
%--------------------

p = inputParser;
p.addParameter('mainstr', '_Scenario',         @(x) validateattributes(x, {'char'}, {}));
p.addParameter('dietstr', '_DC',               @(x) validateattributes(x, {'char'}, {}));
p.addParameter('detfstr', '_DetFate',          @(x) validateattributes(x, {'char'}, {}));
p.addParameter('growstr', '_GrowthParameters', @(x) validateattributes(x, {'char'}, {}));

p.parse(varargin{:});

Opt = p.Results;

%--------------------
% Read files
%--------------------

mainfile = [basename Opt.mainstr '.txt'];
dietfile = [basename Opt.dietstr '.txt'];
detffile = [basename Opt.detfstr '.txt'];
growfile = [basename Opt.growstr '.txt'];

if ~exist(mainfile, 'file')
    error('Main parameter file (%s) not found', mainfile);
end
if ~exist(dietfile, 'file')
    error('Diet file (%s) not found', dietfile);
end
if ~exist(detffile, 'file')
    error('Detritus fate file (%s) not found', detffile);
end


% Main file

fid = fopen(mainfile);
tmp = textscan(fid, '%f', 3);

nrow    = tmp{1}(1);
ncol    = tmp{1}(2);
nstanza = tmp{1}(3);

maindata = textscan(fid, repmat('%f',1,ncol), 'delimiter', '\t');
maindata = cat(2, maindata{:});

fclose(fid);

if size(maindata,1) ~= nrow
    error('#rows value does not match size of data matrix in %s', mainfile);
end

% Diet file

fid = fopen(dietfile);
tmp = textscan(fid, '%f', 2);

nrow = tmp{1}(1);
ncol = tmp{1}(2);

dietdata = textscan(fid, repmat('%f',1,ncol), 'delimiter', '\t');
dietdata = cat(2, dietdata{:});

fclose(fid);

if size(dietdata,1) ~= nrow
    error('#rows value does not match size of data matrix in %s', dietfile);
end

% Detritus fate file

fid = fopen(detffile);
tmp = textscan(fid, '%f', 3);

nrow  = tmp{1}(1);
ncol  = tmp{1}(2);
ndet  = tmp{1}(3);

detfdata = textscan(fid, repmat('%f',1,ncol), 'delimiter', '\t');
detfdata = cat(2, detfdata{:});

if size(detfdata,1) ~= nrow
    error('#rows value does not match size of data matrix in %s', detffile);
end


fclose(fid);

% Growth params file

if nstanza > 0
    if ~exist(growfile, 'file')
        error('Growth parameter file (%s) not found', growfile);
    end

    fid = fopen(growfile);
    tmp = textscan(fid, '%f', 2);

    nrow = tmp{1}(1);
    ncol = tmp{1}(2);

    growdata = textscan(fid, repmat('%f',1,ncol), 'delimiter', '\t');
    growdata = cat(2, growdata{:});

    fclose(fid);

    if size(growdata,1) ~= nrow
        error('#rows value does not match size of data matrix in %s', growfile);
    end

end


%--------------------
% Build ecopathmodel
%--------------------

ngroup = size(maindata,1);
nlive  = ngroup - ndet;
ngear  = 1;

pp = maindata(:,9);
pp(pp == 0) = -1; % Flip consumers and detritus... EwE-F uses reversed convention
pp(pp == 2) = 0;
pp(pp == -1) = 2;

if nstanza > 0
    stzname = strtrim(cellstr(num2str((1:nstanza)', 'stanza%d')));
    EM = ecopathmodel(ngroup, nlive, ngear, 'pp', pp, 'stanzas', stzname);
else
    EM = ecopathmodel(ngroup, nlive, ngear, 'pp', pp);
end

maindata(maindata == -999) = NaN;

% Fill in basic data

EM.groupdata.b        = maindata(:,1);
EM.groupdata.pb       = maindata(:,2);
EM.groupdata.qb       = maindata(:,3);
EM.groupdata.ee       = maindata(:,4);
EM.groupdata.ge       = maindata(:,5);
EM.groupdata.gs       = maindata(:,6);
EM.groupdata.dtImp    = maindata(:,7);
EM.groupdata.ba       = zeros(ngroup,1);
EM.groupdata.immig    = zeros(ngroup,1);
EM.groupdata.emig     = zeros(ngroup,1);
EM.groupdata.stanza   = maindata(:,11);
EM.groupdata.ageStart = maindata(:,12);
EM.groupdata.import(pp == 0) = dietdata(end,:)';

% Diet

dc = zeros(ngroup);
dc(:,pp==0) = dietdata(1:end-1,:);
EM.dc(:,:) = num2cell(dc);

% Flow to detritus

EM.df(:,:) = num2cell(detfdata(:,1:end-1));

% Fisheries (assume discards = 0)

EM.landing(:,:) = num2cell(maindata(:,8));

% Stanza data

EM.stanzadata.stanzaID = (1:nstanza)';
EM.stanzadata.RecPower = growdata(:,2);
EM.stanzadata.BABsplit = growdata(:,3);
EM.stanzadata.WmatWinf = growdata(:,4);

iss = ~isnan(EM.groupdata.stanza);
EM.groupdata.vbK(iss) = growdata(EM.groupdata.stanza(iss),1);







