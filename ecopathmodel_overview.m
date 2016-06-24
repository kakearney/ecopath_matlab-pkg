%% An overview of the ecopathmodel class
%
% The ecopathmodel class offers a Matlab-based implementation of the
% popular Ecopath mass-balance algorithm.  For information on the Ecopath
% concept, I refer you to the official Ecopath with Ecosim website:
% http://ecopath.org/, as well as to the following journal articles:
%
% Christensen, V. & Pauly, D. ECOPATH II--a software for balancing
% steady-state ecosystem models and calculating network characteristics.
% Ecological Modelling 61, 169-185 (1992).  
% 
% Christensen, V. & Walters, C. J. Ecopath with Ecosim: methods,
% capabilities and limitations. Ecological Modelling 172, 109-139 (2004).
%
% This overview assumes that you are already familiar with the Ecopath
% concept, and simply focuses on the use of this particular
% implementation.  Please note that this code is intended to replicate only
% the Ecopath algorithm, not Ecosim, Ecospace, or any of the other
% Ecopath-derived functionalities within the full EwE software.  
%
%% The |ecopathmodel| object
%
% The foundation of this package is the ecopathmodel object, which holds
% all the input data related to a particular modeled ecosystem.  This
% includes the  definition of functional groups and fishing fleets and the
% connections between them, as well as the many parameters associated with
% each state variable and group-to-group flow process.  
%
% Currently, there are three ways to create one of these objects:
%
% # Load data from an EwE6 database file
% # Load data from a set of Rpath-formatted .csv files
% # Build the model manually
%
% If you have a prexisting model, the first two methods are preferable to
% manual transcription, since in my experience it is very easy to
% accidentally mis-transcribe a value, or transcribe it using less
% precision than in the original model.  Building a model manually has the
% advantage of providing a clear "paper trail" for all your model
% parameters, but you have to be careful that you fill in all the necessary
% values.  While I provide a few basic checks of input, I'm not quite as
% good as EwE6 is at  required vs. optional given the architecture of your
% specific ecosystem.
%
% Personally, I prefer a combination of the two techniques, where possible.
% Use the EwE6 software during the initial data-gathering step of building
% an Ecopath model.  That will allow you to take advantage of all the
% graphical utilities that warn you when parameters are incorrect.  Once
% the model data is acceptable, and the Ecopath Parameterization
% calculations in EwE6 return a balanced (or at least nearly-balanced)
% result, move over to this tool for project-specific parameter
% adjustments.  This process allows you to preserve one copy of the base
% model while also keeping project-specific details clearly documented.

%% Importing Ecopath with Ecosim data
%
% The current version of Ecopath with Ecosim stores data in
% specially-formatted Microsoft Access database files, with the .EwEmdb
% extension.  The |mdb2ecopathmodel| function will read data from one of
% these files into an |ecopathmodel| object, with two caveats:
%
% # You need to first install the |mdbtools| set of command-line utilities.
% They're available for download here: https://github.com/brianb/mdbtools.
% Install instruction for *nix systems are included on GitHub; for Mac, I
% recommend installing via either MacPorts or Homebrew.  Windows users (and
% I know that's most of the Ecopath-using community): I've never compiled
% this tool on a Windows machine, though I know it can be done.  I don't
% work with any Windows machines myself, so this is a big shortcoming in
% this code base at the moment. If any of you have compiled mdbtools on a
% Windows system, or have a better way of reading .mdb files into Matlab on
% a Windows machine (using the Database Toolbox, maybe?), please contact
% me!  I'd love to make this utility easier to use on Windows.
% # I expect all files to use the EwE6 format.  If you have older files,
% you'll first need to import them into EwE6 and let its conversion tool
% convert to the newer format.
%
% This example reads in one of the example ecosystems that ships with the
% EwE6 software:
%

epfolder = '~/Documents/Research/Working/EcopathModels/';
Gen37 = mdb2ecopathmodel(fullfile(epfolder, 'Generic_37.EwEmdb'));

%%
% When importing, you'll usually encounter a few warnings like those seen
% above.  EwE6 uses a few NaN-placeholders in its data files, which are
% then converted to 0s when it performs the Ecopath calculation.  My code
% replaces those placeholders on reading, and lets you know.  It also
% alters names, if necessary.  If you're not happy with the "translation", 
% you can  alter the |name|, |fleet|, or |stanza| properties of the
% |ecopathmodel| object, and these changes will propagate to the rest of
% the tables.  For example, I don't like those trailing underscores on a
% couple groups:   

Gen37.name{11} = 'Pelagics_Small_Carniv';
Gen37.name{12} = 'Pelagics_Small_Herbiv';

Gen37.groupdata

%% Importing Rpath data
%
% Rpath is an R-based implementation of Ecopath with Ecosim, written by Sea
% Lucey and Kerim Aydin.  It is available for download on GitHub:
% https://github.com/slucey/RpathDev.  The Rpath package includes two
% functions, read.rpath.params and write.rpath.params, to import and export
% data from comma-delimited text files.
%
% This code includes a function to read data from .csv files that match the
% format used by those two functions.  The following example reads in the
% model described in the Rpath vignette (REco.params):

rfolder = '~/Documents/Research/Working/Rpath/tests/';
REco = rpath2ecopathmodel(fullfile(rfolder, 'REco'));

%%
% Again this function will typically issue several warinngs related to the
% differing ways this code and Rpath use NaNs vs 0s as placeholders for
% certain parameters.  Rpath also assigns pedigree values to all groups,
% even non-leading stanza groups; my code doesn't allow that so those
% values are stripped out of the pedigree table.


%% Building an |ecopathmodel| object manually
%
% As I mentioned above, I don't really recommend that you start a model
% from scratch using just this tool.  My focus when developing this code
% was to increase the flexibility of the Ecopath algorithm, not to
% replicate the EwE6 GUI capabilities (including its many data validation
% checks).  
%
% However, Ecopath models are very often published in peer-reviewed
% journals or technical reports without accompanying data files.  In these
% cases, you many need to manually transcribe the data from various printed
% tables in order to carry out additional calculations.
%
% For this example, I'm going to use the Eastern Pacific Subarctic Gyre
% ecosystem model, which is the one I happened to be working with during
% the development of this class, and which required manual transcription.
% The details of that model were published in
%
% Aydin KY, McFarlane GA, King JR, Megrey BA (2003) The BASS/MODEL report
% on trophic models of the Subarctic Pacific Basin ecosystems. PICES Sci
% Rep 25  

%%
% Start by creating an empty |ecopathmodel| object.  For this, you need 3
% parameters: 
%
% * number of total groups
% * number of live groups
% * number of fishing gears/fleets.
%
% While not required, it's highly recommended that you also add the names
% of all groups, fleets, and stanzas, since that data is used to set up and
% label all the table columns and rows, and the type of each group (i.e.
% whether consumer, producer, or mixotroph), since that is used to validate
% other parameters' values as they're added.  

names = {...
'Sperm whales'
'Toothed whales'
'Fin whales'
'Sei whales'
'Northern fur seals'
'Elephant seals'
'Dall''s porpoises'
'Pacific white-sided dolphins'
'Northern right whale dolphins'
'Albatross'
'Shearwaters'
'Storm Petrels'
'Kittiwakes'
'Fulmars'
'Puffins'
'Skuas'
'Jaegers'
'Sharks'
'Large gonatid squid'
'Boreal clubhook squid'
'Neon flying squid'
'Sockeye salmon'
'Chum salmon'
'Pink salmon'
'Coho salmon'
'Chinook salmon'
'Steelhead'
'Pomfret'
'Saury'
'Pelagic forage fish'
'Micronektonic squid'
'Mesopelagic fish'
'Large jellyfish'
'Ctenophores'
'Salps'
'Chaetognaths'
'Sergestid shrimp'
'Misc predatory zooplankton'
'Amphipods'
'Pteropods'
'Euphausiids'
'Copepods'
'Microzooplankton'
'Bacteria'
'Large phytoplankton'
'Small phytoplankton'
'DNH3'
'POM'};

ngroup = length(names);
ngear = 1;          % Ecopath requires at least 1, even if it catches nothing
nlive = ngroup - 2; % DNH3 and POM are detrital
isprod = ismember(names, {'Large phytoplankton', 'Small phytoplankton'});
isdet = ismember(names, {'DNH3', 'POM'});

pp = zeros(ngroup,1); % 0 = consumer
pp(isprod) = 1;       % 1 = producer
pp(isdet) = 2;        % 2 = detritus

% Names of groups, fleets, and stanzas must meet Matlab's variable name
% restrictions, since they will be used as table row/column labels.  To
% meet this requirement, here I capitalize all words and then strip out
% spaces and special characters. 

names = regexprep(names,'(\<[a-z])','${upper($1)}');
names = regexprep(names, '[\s-\.'']', '');

% Now create empty ecopathmodel

Esa = ecopathmodel(ngroup, nlive, ngear, 'groups', names, 'pp', pp)

%%
% The ecoathmodel object properties include several table arrays.  You can
% refer to the ecopathmodel property descriptions (|help ecopathmodel|) to
% see what parameters are stored in each table.  The variable names are all
% based on those used in EwE6, so they should be familiar to most Ecopath
% users.    

%%
% Now it's time to start adding the data.  We'll start with the groupdata
% table, which holds all the group-related parameters.

Esa.groupdata

%%
% The first several columns of the groupdata table correspond to the "Basic
% input" panel in EwE6.  These, along with diet fractions, are the
% parameters most likely to be published in any Ecopath-related study.
% In this case, the values come straight from Table B6 in the Aydin et al.,
% 2003 report.
%
% You'll notice a few warning messages printed to the screen as I add the
% data into the appropriate tables; those are letting me know that I tried
% to add invalid values to certain locations.  In this case, the
% discrepancy is that the printed table left some values blank where
% internally they're actually supposed to be 0; those incorrect NaNs are
% replaced with the appropriate 0s by the data validator.

% Columns are trophic level (TL), biomass (B, t/km2), production/biomass
% (P/B, 1/year), consumption/biomass (Q/B, 1/year), ecotrophic efficiency
% (EE, proportion), growth efficiency (PC, proportion) biomass accumulation
% (BA t/km2/year), unassimilated respiration (UnAss, proportion) and the
% proportion of detritus flowing to NH3 and POM, respectively

tableB6 = [...
5.4 0.000929    0.0596  6.61    0       0.00902 0 0.2 0.5 0.5
5.2 0.000028    0.0252  11.16   0       0.00226 0 0.2 0.5 0.5
4.1 0.027883    0.02    4.56    0.12912 0.00439 0 0.2 0.5 0.5
4.1 0.005902    0.02    6.15    0.1358  0.00325 0 0.2 0.5 0.5
5.2 0.000246    0.235   39.03   0.01083 0.00602 0 0.2 0.5 0.5
5.2 0.00043     0.368   11.08   0.00692 0.03321 0 0.2 0.5 0.5
5.2 0.00598636  0.1     27.47   0.02546 0.00364 0 0.2 0.5 0.5
5.2 0.00396248  0.14    25.83   0.01819 0.00542 0 0.2 0.5 0.5
5.3 0.00389728  0.16    24.14   0.01592 0.00663 0 0.2 0.5 0.5
5.9 0.00004     0.05    81.59   0.05043 0.00061 0 0.2 0.5 0.5
4.7 0.0004      0.1     100.13  0.02547 0.001   0 0.2 0.5 0.5
4.6 0.000056    0.1     152.08  0.02546 0.00066 0 0.2 0.5 0.5
4.6 0.000052    0.1     123     0.02549 0.00081 0 0.2 0.5 0.5
4.9 0.000074    0.1     100.26  0.02557 0.001   0 0.2 0.5 0.5
4.7 0.000058    0.1     104.33  0.02535 0.00096 0 0.2 0.5 0.5
4.8 0.000054    0.075   96.6    0.0338  0.00078 0 0.2 0.5 0.5
4.8 0.000038    0.075   96.6    0.03388 0.00078 0 0.2 0.5 0.5
5.4 0.05        0.2     10.95   0       0.01826 0 0.2 0.5 0.5
4.2 0.03        2.555   7.3     0.19453 0.35    0 0.2 0.5 0.5
4.9 0.012       2.555   7.3     0.19453 0.35    0 0.2 0.5 0.5
5.3 0.45        2.555   6.205   0.91095 0.41176 0 0.2 0.5 0.5
4.3 0.08965573  1.27    10.13   0.32249 0.12537 0 0.2 0.5 0.5
3.7 0.05413587  1.93    14.51   0.21221 0.13301 0 0.2 0.5 0.5
4.2 0.02326662  3.37    18.49   0.12153 0.18226 0 0.2 0.5 0.5
4.9 0.00445349  2.47    16.55   0.16581 0.14924 0 0.2 0.5 0.5
4.9 0.00930315  0.8     5.33333 0.51195 0.15    0 0.2 0.5 0.5
4.9 0.0093      0.8     5.33333 0.51212 0.15    0 0.2 0.5 0.5
4.8 0.21        0.75    3.75    0.54697 0.2     0 0.2 0.5 0.5
3.8 0.45        1.6     7.9     0.5545  0.20253 0 0.2 0.5 0.5
3.9 0.92156     1.5     5       0.9     0.3     0 0.2 0.5 0.5
3.9 0.87135     3       15      0.9     0.2     0 0.2 0.5 0.5
3.9 4.5         0.9     3       0.16002 0.3     0 0.2 0.5 0.5
3.6 4           3       10      0       0.3     0 0.2 0.5 0.5
2.7 9.1         4       110     0.05269 0.03636 0 0.2 0.5 0.5
2.7 8           9       30      0.02371 0.3     0 0.2 0.5 0.5
3.5 6.6         2.555   12.045  0.27638 0.21212 0 0.2 0.5 0.5
3.5 5           2.555   12.045  0.18813 0.21212 0 0.2 0.5 0.5
3.5 5.0688      2.555   12.045  0.18871 0.21212 0 0.2 0.5 0.5
3.1 10.1376     2.555   12.045  0.64429 0.21212 0 0.2 0.5 0.5
3.1 10.1376     2.555   12.045  0.53491 0.21212 0 0.2 0.5 0.5
3.1 25.344      2.555   12.045  0.53934 0.21212 0 0.2 0.5 0.5
2.4 34.848      23.725  112.42  0.88106 0.21104 0 0.2 0.5 0.5
2.3 35          48.91   233.235 0.99619 0.2097  0 0.2 0.5 0.5
2   122.9031    18.45   25      0.9     0.738   0 0.2 0.5 0.5
1   69.7        42.34   NaN     0.67337 NaN     0 NaN 0.5 0.5
1   76          129.575 NaN     0.77256 NaN     0 NaN 0.5 0.5
1   NaN         NaN     NaN     0.42757 NaN     0 NaN 0   0
1   NaN         NaN     NaN     0.42757 NaN     0 NaN 0   0];  

% Plug these values into the appropriate groupdata columns:

Esa.groupdata.b  = tableB6(:,2);
Esa.groupdata.pb = tableB6(:,3);
Esa.groupdata.qb = tableB6(:,4);
Esa.groupdata.ee = tableB6(:,5);
Esa.groupdata.ge = tableB6(:,6);
Esa.groupdata.ba = tableB6(:,7);
Esa.groupdata.gs = tableB6(:,8);

% The detrital flows (last two columns of the published table) go into the
% df table: 

Esa.df.DNH3 = tableB6(:,9);
Esa.df.POM  = tableB6(:,10);

% Certain values in that table were marked in gray, indicating that they
% were calculated by the Ecopath mass-balance calculation, not provided as
% input.  We'll get rid of those values for now:

nob  = [30 31 44];
noqb = 26:28;
noee = [1:29 32:43 45:48];
noge = [1:25 29:48];

Esa.groupdata.b(nob)   = NaN;
Esa.groupdata.bh(nob)  = NaN;
Esa.groupdata.qb(noqb) = NaN;
Esa.groupdata.ee(noee) = NaN;
Esa.groupdata.ge(noge) = NaN;

%%
% Next up, diet fractions.  Here's one of my personal pet peeves:
% publishing sparse matrices (like Ecopath diet fraction matrices) as  
% tables of mostly empty space spread across multiple printed pages.  Can
% we all agree to stop doing this?  Pretty please?     
%
% Anyway, diet data was copied and pasted from Table B8 (all 5 pages of it)
% of the report into a csv file, which I've pasted here to minimize the
% number of supporting files needed for this example.      

dietcsv = {...
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.22725,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.04811,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.002,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00351,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.04879,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.03229,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.03176,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00032,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00326,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00046,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00042,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00061,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00047,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00044,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',0.00031,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.02287,0.00305,0.001524,0.001524,0.00915,0.0122,0.01372,0.00762,0.01524,0.04573,,,,,,,,0.01782,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.00915,0.00122,0.000609,0.000609,0.00366,0.00488,0.00549,0.00305,0.0061,0.01829,,,,,,,,0.00713,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.34299,0.04573,0.022866,0.022866,0.1372,0.18293,0.20579,0.11433,0.22866,0.68598,,,,,,,,0.26724,,,0.295,,,,,,,,,,,,,,,,,,,,,,,'
'0.00208,0.0249,0.001245,0.001245,0.02988,0.01132,0.01494,0.03984,0.00498,,,,,,,,,0.05324,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.00125,0.01504,0.000752,0.000752,0.01804,0.00683,0.00902,0.02406,0.00301,,,,,,,,,0.03215,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.00054,0.00646,0.000323,0.000323,0.00775,0.00294,0.00388,0.01034,0.00129,,,,,,,,,0.01382,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.0001,0.00124,6.19E-05,6.19E-05,0.00148,0.00056,0.00074,0.00198,0.00025,,,,,,,,,0.00264,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.00022,0.00258,0.000129,0.000129,0.0031,0.00117,0.00155,0.00413,0.00052,,,,,,,,,0.00552,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.00022,0.00258,0.000129,0.000129,0.0031,0.00117,0.00155,0.00413,0.00052,,,,,,,,,0.00552,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.00486,0.05833,0.002917,0.002917,0.06999,0.056,0.035,0.09332,0.01167,,,,,,,,,0.12471,,,,,,,,,,,,,,,,,,,,,,,,,,'
'0.01042,0.12498,0.00625,0.00625,0.14998,0.12,0.07499,0.19997,0.025,0.1,0.275,0.05,0.4,,0.4,0.5,0.5,0.26724,,,0.058,,,,,,,0.04,,,,,,,,,,,,,,,,'
'0.126,0.264,0.076,0.076,0.267,0.2,0.208,0.172,0.053,0.1,0.275,0.05,0.4,0.04,0.4,0.5,0.5,0.103,0.01,0.01,0.319,0.10982,0.00802,0.068096,0.367206,0.367206,0.367206,,,,,,,,,,,,,,,,,'
'0.375,0.05,0.025,0.025,0.15,0.4,0.225,0.225,0.25,0.05,0.3,0.6,,0.96,0.1,,,0.1,0.33,0.99,0.223,0.07968,0.03929,0.034823,0.205691,0.205691,0.205691,0.75,,,0.05,,,,,,,,,,,,,'
'0.10416,,0.062499,0.062499,0.15,,0.2,0.1,0.4,,,,,,,,,,,,0.105,0.10982,0.00802,0.068096,0.367206,0.367206,0.367206,0.08,,,,,,,,,,,,,,,,'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
',,,,,,,,,,,,,,,,,,,,,0.01609,0.20285,0.003826,7.31E-05,7.31E-05,7.31E-05,,,,,,0.04356,,,,,,,,,,,'
',,,,,,,,,,,,,,,,,,,,,0.01609,0.20285,0.003826,7.31E-05,7.31E-05,7.31E-05,,,,,,0.03829,,,,,,,,,,,'
',,0.054357,0.05435,,,,,,,,,,,,,,,0.04484,,,,0.00043,,,,,0.01,0.05298,0.06795,0.06455,0.15,0.03159,,,,,,,,,,,'
',,0.041179,0.041179,,,,,,,,,,,,,,,0.03397,,,,,,,,,0.01,0.04014,0.05147,0.0489,0.03,0.02393,,,,,,,,,,,'
',,0.041746,0.041746,,,,,,,,,,,,,,,0.03444,,,0.00171,0.01474,8.64E-06,,,,0.01,0.04069,0.05218,0.04957,0.03,0.02426,,,,,,,,,,,'
',,0.083492,0.083492,,,,,,,0.0189,0.0378,0.0252,,0.0126,,,,0.06888,,,0.29328,0.07363,0.32162,0.008047,0.008047,0.008047,0.03,0.08138,0.10437,0.09915,0.24,0.04852,,,0.04444,0.04444,0.04444,,,,,,'
',,0.083492,0.083492,,,,,,,0.0189,0.0378,0.0252,,0.0126,,,,0.06888,,,0.23976,0.06719,0.440933,0.04623,0.04623,0.04623,0.01,0.08138,0.10437,0.09915,0.031,0.04852,,,0.04444,0.04444,0.04444,,,,,,'
',,0.20873,0.20873,,,,,,,0.04724,0.09449,0.06299,,0.0315,,,,0.1722,,,0.10058,0.08491,0.038951,0.002736,0.002736,0.002736,0.05,0.20344,0.26091,0.24787,0.171,0.12131,,,0.11111,0.11111,0.11111,,,,,,'
',,0.287004,0.287004,,,,,,,0.06496,0.12992,0.08661,,0.04331,,,,0.23678,,,0.03318,0.29806,0.019816,0.002736,0.002736,0.002736,0.01,0.5,0.35875,0.34082,0.348,0.62,0.25,0.25,0.8,0.8,0.8,0.4,0.4,0.4,,,'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,0.25,0.25,,,,0.4,0.4,0.4,0.3,,'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,0.25,0.25,,,,,,,,0.25,'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,0.25,0.25,,,,0.2,0.2,0.2,0.4,,'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,0.3,0.75,'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,0.5'
',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,0.5'
};

dcfile = [tempname '.csv'];
fid = fopen(dcfile, 'w+');
fprintf(fid, '%s\n', dietcsv{:});
fclose(fid);

dc = csvread(dcfile);

delete(dcfile);

%%
% One of the reasons pdf tables are a terrible way to share data: rounding
% error.  In this case, some values must have been rounded for display;
% diet data (plus import) should always sum to 1 for each predator, and in
% this case it doesn't.  

sum(dc,1)

%%
% We'll normalize to get a proper diet fraction table.

dc = bsxfun(@rdivide, dc, sum(dc,1));

Esa.dc(:,:) = num2cell(dc);


%%
% There are a few remaining variables that weren't really discussed in the
% Aydin et al. report, but are necessary to replicate their Ecopath
% calculations.

% Detritus import: assume to be 0 for both detritus groups

Esa.groupdata.dtImp(isdet) = zeros(ngroup-nlive,1);

% Immigration/emigration: Assume 0

Esa.groupdata.immig = zeros(ngroup,1);
Esa.groupdata.emig  = zeros(ngroup,1);

% Fisheries landings and discards: Assume no fisheries loss in this model.

Esa.landing(:,:) = num2cell(zeros(ngroup,ngear));
Esa.discard(:,:) = num2cell(zeros(ngroup,ngear));

% Discard fate: Values don't really matter, since there are no discards to
% redirect, but as a placeholder, split this between the two detrital
% groups.

Esa.discardFate(:,:) = {0.5, 0.5};

% Finally, the detritus groups are currently listed with no biomass or
% import.  The Aydin et al., 2003 work used a precursor to the Rpath model
% for its calculations, rather than the more familiar EwE platform. One
% difference between Rpath and EwE is that Rpath calculates detrital pool
% biomass based on an assumed overturning rate.  This Matlab package also
% allows that possibility, via the |detpb| column.  However, the rate used
% in this particular study wasn't published, so I took the easy way out and
% simply substituted some approximate values (for Ecopath balance purposes,
% this number doesn't affect the calculations at all).

Esa.groupdata.b(isdet) = ones(ngroup-nlive,1)*50;

%% 
% At this point, we have a fully-populated ecopathmodel object


%% Filling in multi-stanza group parameters
% 
% Multi-stanza groups in an Ecopath model represent different life stages
% of a single functional group.  In time-dynamic models such as Ecosim,
% growth leads to a flux of biomass from younger groups to older groups.
% Although this tool does not explicitly include Ecosim-like calculations,
% I do strive to maintain consistency with EwE6 and Rpath, which means that
% multistanza groups must match the stable age distibution requirements.
%
% The primary |ecopathmodel| method for calculating stanza-related
% parameters is the |calcstanza| method, which fills in the biomass (b), 
% comsumption rates (qb), and biomass accumulation rates (ba) for all
% groups that are part of a multistanza set.
%
% As an example, let's look at the REco model that we read in above.  This
% model includes 4 multi-stanza sets, with two groups apeice:

REco.stanzadata
REco.groupdata

%%
% As you can see, the B and QB values for all the juvenile groups are
% currently set to NaN.  Running the |calcstanza| function will fill these
% in:

REco = REco.calcstanza;
REco.groupdata

%%
% You can look at the stable growth curve plots by using the |'plot'|
% option:

REco.calcstanza('plot', true);

%%
% The |setstanzas| method does almost the same thing as the |calcstanza|
% method.  However, it provides a few additional checks before filling in
% the B and QB values.  This is advantageous if you import data from an
% EwE6 file, which has the non-leading stanza parameters already filled in.
%  For example, let's look at another stock example from the EwE6 software,
%  the Tampa Bay model:

Tb = mdb2ecopathmodel(fullfile(epfolder, 'Tampa_Bay.EwEmdb'));
Tb.groupdata
Tb.stanzadata

%%
% As you can see, the non-leading group data is already present, so we
% don't really need to recalculate it. But what would happen if we did? 

Tb2 = Tb.calcstanza;

bvals  = [Tb.groupdata.b  Tb2.groupdata.b  Tb.groupdata.b  - Tb2.groupdata.b]
qbvals = [Tb.groupdata.qb Tb2.groupdata.qb Tb.groupdata.qb - Tb2.groupdata.qb]

%%
% The values aren't exactly the same.  The differences arise because of
% some discrepancies in the way the growth curve calculations are done in
% my code vs. in EwE6...specifically, how the tail end (age of last 10% of
% biomass to infinity) is handled (Rpath handles things the same as I do).
% The differences will be more noticeable in long-lived stanza sets,
% especially in the QB variable.
%
% The |setstanzas| method keeps an eye out for little differences between
% parameters that are already filled in and ones it tries to calculate.  If
% it finds a difference that is less than 0.5% of the value, it won't alter
% the already-filled in data.  If it finds a larger difference, it will
% change it, but it will issue a warning (this is likely an indication that
% the original data source used a different calculation for the stable
% growth curve).

Tb2 = Tb.calcstanza;
Tb3 = Tb.setstanzas;
isequaln(Tb.groupdata.b, Tb2.groupdata.b)
isequaln(Tb.groupdata.b, Tb3.groupdata.b)









