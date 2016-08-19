classdef ecopathmodel
%ECOPATHMODEL A Matlab-based version of the Ecopath food web model
%   
% The ecopathmodel class forms the basis for my Matlab-based implementation
% of Ecopath, adapted from the Ecopath with Ecosim (EwE) model
% (www.ecopath.org). 
%
% The EwE model is a ecosystem model used primarily in the fisheries
% modeling community.  In this community, the term "Ecopath model" is used
% in two contexts: 
% - To refer to the algorithm that estimates a snapshot of an ecosystem
%   (including the amount of biomass in various ecosystem components, and
%   the fluxes between these components) 
% - To refer to a set of input data for a particular ecosystem, intended to
%   be used with the Ecopath (and/or Ecosim, Ecospace, etc.) algorithm(s)
%
% In this code, the ecopathmodel object represents the latter of these: a
% storage-container for input data related to the Ecopath algorithm.  
%
% The units used in the property descriptions are defined in terms of three
% dimensions: M = mass, A = area (or volume), and T = time. In the original
% software, the default is M = metric tons wet weight, A = km^2, and T =
% year.  You can use whatever units work best for your own purposes, as
% long as they remain consistent across all variables.      
% 
% Constructor method syntax:
%
%   obj = ecopathmodel(ngroup, nlive, ngear, p1, v1, ...)
%
% Please type "help ecopathmodel.ecopathmodel" for full syntax details
% 
% ecopathmodel Properties:
%
%   ngroup:     scalar, Number of groups (living and detrital) in the model
%
%   nlive:      scalar, Number of living groups (non-detrital) in model
%
%   ngear:      scalar, Number of fishing gears/fleets in model
%
%   name:       ngroup x 1 cell array of strings, names corresponding to
%               each group in the model. 
%
%   fleet:      ngear x 1 cell array of strings, names corresponding to
%               each fishing gear in the model 
%
%   groupdata:  table of parameters related to functional groups.  Rows
%               correspond to functional groups; columns as follows: 
%
%               b:          biomass (M A^-1).  A NaN indicates a value to
%                           be estimated by the ecopath algorithm.           
%
%               pb:         production/biomass ratio (T^-1). A NaN
%                           indicates a value to be estimated by the
%                           ecopath algorithm.  
%
%               qb:         consumption/biomass ratio (T^-1). A NaN
%                           indicates a value to be estimated by the
%                           ecopath algorithm.  
%
%               ee:         ecotrophic efficiency (unitless). A NaN
%                           indicates a value to be estimated by the
%                           ecopath algorithm.  
%
%               ge:         growth efficiency (unitless). A NaN indicates a
%                           value to be estimated by the ecopath algorithm. 
%
%               gs:         fraction of consumed food that is not
%                           assimilated (unitless)
%
%               dtImp:      detritus import (should be zero for all
%                           non-detrital groups) (MA^-1 T^-1) 
%
%               bh:         habitat biomass, i.e. biomass per unit
%                           habitable area (M A^-1).  This variable is
%                           designed as a shortcut if all your critters are
%                           clustered in only a small portion of the
%                           habitat area, such that bh = b/areafrac.  For
%                           each group, either b or bh can be defined;
%                           leave as NaN to calculate based on the other.
%
%               pp:         fraction of diet consisting of primary
%                           production. pp = 0 indicates a strict consumer,
%                           pp = 1 is a primary producer, pp = 2 indicates
%                           a detrital group, 0<pp<1 indicates a mixotroph.
%
%               areafrac:   fraction of habitat area occupied by each group
%                           (unitless, 0-1)  
%
%               ba:         biomass accumulation (M A^-1 T^-1)
%
%               baRate:     biomass accumulation per unit biomass (T^-1).
%                           For each group, a value may be entered for
%                           either ba or baRate; leave as NaN to calculate
%                           based on the other.
%
%               immig:      immigration into area (M A^-1 T^-1)  
%
%               emig:       emigration out of area (M A^-1 T^-1) 
%
%               emigRate:   emigration per unit biomass (T^-1).  For each
%                           group, a value may be entered for either emig
%                           or emigRate; leave as NaN to calculate based on
%                           the other.    
%
%               stanza:     index of stanza set to which each group
%                           belongs.  A value of 0 indicates that the group
%                           is not part of a multi-stanza set.
%
%               ageStart:   Age (in months) at which this group begins.
%                           Will be 0 for all non-stanza groups.
%
%               vbK:        K parameter from von Bertalanffy growth curve,
%                           applicable to groups that are part of a
%                           multi-stanza set.  Should be same value for all
%                           age groups within a single stanza set.  NaN for
%                           all non-stanza groups.
%
%               detpb:      detrital production rate (T^-1). If not NaN,
%                           (for detrital groups only), indicates that the
%                           biomass for this group should be estimated
%                           assuming a turnover rate of the detrital pool
%                           equal to this value. 
%
%               import:     fraction of group's diet that comes from
%                           outside the system (unitless, 0-1)
%
%   dc:         table of diet fraction data. Rows correspond to prey
%               groups, columns to predators.  The sum of each column (plus
%               import) should be 1.
%
%   landing:    table of landing data (M A^-1 T^-1).  Fisheries landings of
%               each group by each gear.  Rows correspond to functional
%               groups, columns to gear types.
%
%   discard:    table of discard data (M A^-1 T^-1).  Fisheries discards of
%               each group by each gear.  Rows correspond to functional
%               groups, columns to gear types.
%
%   df:         table of detritus fate (unitless, 0-1).  Values represent
%               fraction of non-predatory losses (unassimilated
%               consumption, non-predatory mortality) that go into each
%               detrital pool. Rows correspond to functional groups,
%               columns to detrital groups.
%
%   discardFate:table of discard fate (unitless, 0-1).  Values represent
%               fraction of fisheries discards that go into each detrital
%               pool. Rows correspond to gear types, columns to detrital
%               groups.   
%
%   stanza:     nstanza x 1 cell array of strings, names corresponding to
%               stanza sets.  Stanza sets are used to indicate that 2 or
%               more functional groups represent different life stages of a
%               single species/classification. (In EwE, these are sometimes
%               referred to as multi-stanza groups.  I find the use of the
%               word "group" a little confusing, since it sometimes refers
%               to the all-ages parent grouping and sometimes to the
%               individual sub-groups, so in this code I use the term
%               "stanza set" for the parent container and "stanza group"
%               for a functional group that is part of a set).
%
%   stanzadata: table holding parameters related to each multi-stanza set.
%               Rows correspond to stanza sets; columns are as follows:
%
%               stanzaID:   index of each stanza set
%
%               BABsplit:   biomass accumulation rate for the stanza set
%                           (T^-1).  
%
%               Btot:       total biomass across all stanza groups in a
%                           set.  Can be used as an alternative reference
%                           (as opposed to leading stanza biomass) when
%                           calculating biomass-per-stanza-group curves.
%
%   pedigree:   table of pedigree values for a model.  These values
%               represent the uncertainty of a value as a fraction of that
%               value. Pedigree values can be applied to any element in the
%               ecopathmodel table properties.  Columns are as
%               follows:
%           
%               property:   name of table (groupdata, dc, landing,
%                           discard, df, discardFate, stanzadata)
%
%               row:        row index
%
%               column:     column index
%
%               pedigree:   pedigree value for the parameter stored in
%                           EM.property(row,column)
% 

% Copyright 2016 Kelly Kearney
% kakearney@gmail.com

    
    properties (SetAccess = immutable)
        %NGROUP Number of groups (living and detrital) in model
        ngroup;
        %NLIVE Number of living groups (non-detrital) in model
        nlive;
        %NGEAR Number of fishing gears/fleets in model
        ngear;
    end
    properties
        %NAME Names corresponding to each group in the model
        name;
        %FLEET Names corresponding to each fishing gear/fleet in the model
        fleet;
        %GROUPDATA Table of group-related parameters
        groupdata;
        %DC Table of diet composition data
        dc;
        %LANDING Table of fisheries landings
        landing;
        %DISCARD Table of fisheries discards
        discard;
        %DF Table of detritus fates
        df;
        %DISCARDFATE Table of discard fates
        discardFate;
        %STANZA Names corresponding to each multi-stanza set
        stanza;
        %STANZADATA Table of multi-stanza-set-related parameters
        stanzadata;
        %PEDIGREE Table of pedigree values applied to parameters
        pedigree;
    end

    methods
        
        %------------------------------------------------------------------
        % Constructor
        function obj = ecopathmodel(ngroup, nlive, ngear, varargin)
        %ECOPATHMODEL Create an ecopathmodel object
        %
        % obj = ecopathmodel(ngroup, nlive, ngear, p1, v1, ...)
        %
        % This function creates an ecopathmodel object, setting up all
        % tables to the proper size but with NaNs as placeholders for data.
        % 
        % Input variables:
        %
        %   ngroup:     number of functional groups to include in the
        %               model.  Must be greater than 1.
        %   
        %   nlive:      number of live functional groups in the model.
        %               Must be between 1 and ngroup-1 (all models must
        %               include at least one living group and at least one
        %               detrital group).   
        %
        %   ngear:      number of fishing fleets/gears.  At least 1 gear
        %               must be specifed (landings and discards for this
        %               gear may be set to 0 later to simulate no
        %               fisheries).   
        %
        % Optional input variables, passed as parameter/value pairs:
        %
        %   groups:     ngroup x 1 cell array of strings, names
        %               corresponding to each group.  If not included,
        %               defaults will be groupX, where X is the group
        %               index.  Note that detrital groups should be listed
        %               last.  Because these will be used as table
        %               row/column headers, these strings should conform to
        %               Matlab's restrictions for a valid variable name.
        %
        %   fleets:     ngear x 1 cell array of strings, names
        %               corresponding to each fishing gear. Because these
        %               will be used as table row/column headers, these
        %               strings should conform to Matlab's restrictions for
        %               a valid variable name.   
        %
        %   stanzas:    nstanza x 1 cell array of strings, where nstanza
        %               refers to the number of multi-stanza sets of
        %               groups.  Functional groups that are part of a
        %               multi-stanza set have certain constraints set on
        %               their B, P/B, and Q/B values. Because these will be
        %               used as table row/column headers, these strings
        %               should conform to Matlab's restrictions for a valid
        %               variable name.
        %
        %   pp:         ngroup x 1 array, functional group types, where 0 =
        %               consumer, 1 = producer, 0 < x < 1 = partial primary
        %               producer, 2 = detritus.  Note that detrital groups
        %               are always listed last, so pp((nlive+1):ngroup)
        %               must be equal to 2.  If not entered, the parameters
        %               are initialized assuming groups 1:nlive are
        %               consumers and (nlive+1):ngroup are detrital. 
            
            % Check that ngroup, nlive, and ngear are postive integers
        
            validateattributes(ngroup, {'numeric'}, {'scalar', 'positive'}, '', 'ngroup');
            validateattributes(nlive, {'numeric'}, {'scalar', 'positive', '<', ngroup}, '', 'nlive');
            validateattributes(ngear, {'numeric'}, {'scalar', 'positive'}, '', 'ngear');
            
            if any(rem([ngroup, nlive, ngear],1))
                error('ngroup, nlive, and ngear must be integer values (of any numeric class)');
            end
            
            % Parse inputs and set defaults if not included
            
            grpdefault = strtrim(cellstr(num2str((1:ngroup)', 'group%d')));
            fltdefault = strtrim(cellstr(num2str((1:ngear)', 'fleet%d')));
            stzdefault = {};
            ppdefault = [zeros(nlive,1); ones(ngroup-nlive,1)*2];

            p = inputParser;
            addParameter(p, 'groups', grpdefault, @(x) iscell(x) && length(x)==ngroup && all(cellfun(@ischar, x)));
            addParameter(p, 'fleets', fltdefault, @(x) iscell(x) && length(x)==ngear && all(cellfun(@ischar, x)));
            addParameter(p, 'stanzas',stzdefault, @(x) iscell(x) && (isempty(x) || all(cellfun(@ischar, x)))); 
            addParameter(p, 'pp', ppdefault, @(x) validateattributes(x, {'numeric'}, {'>=', 0, '<=', 2}));
            parse(p, varargin{:});
            
            obj.ngroup = ngroup;
            obj.nlive = nlive;
            obj.ngear = ngear;
            
            % Check that there is at least one detrital group
            
            if obj.nlive > obj.ngroup-1
                error('ecopathmodel must include at least one detrital group (nlive < ngroup-1)')
            end
            
            % Make sure detrital groups are listed last in the pp array
            
            if any(p.Results.pp > 1 & p.Results.pp < 2)
                error('pp values must be in range 0<=pp<=1 or pp == 2');
            end
            if ~all(p.Results.pp((nlive+1):ngroup) == 2)
                error('pp values of last ngroup-nlive groups must be 2 (i.e. detrital groups last)');
            end
            
            % Set default data values (the input checks will add the
            % appropriate placeholders based on group type) 
            
            gvars = {'b', 'pb', 'qb', 'ee', 'ge', 'gs', 'dtImp', 'bh', ...
                    'pp', 'areafrac', 'ba', 'baRate', 'immig', 'emig', ...
                    'emigRate', 'stanza', 'ageStart', 'vbK', 'detpb', ...
                    'import'};
            nv = length(gvars);
            
            ndet = ngroup - nlive;
            
            obj.name = p.Results.groups;
            obj.fleet = p.Results.fleets;
            obj.stanza = p.Results.stanzas;
            
            nstan = length(obj.stanza);
            
            gdtmp = nan(ngroup,nv);
            ispp = strcmp(gvars, 'pp');
            gdtmp(:,ispp) = p.Results.pp;
            gdtmp(:,end) = 0; % import
            
            W = warning('off', 'Ecopathmodel:groupdataValidation');
            
            obj.groupdata = array2table(gdtmp, ...
                'VariableNames', gvars, 'RowNames', obj.name);
            
            obj.dc = array2table(zeros(ngroup), ...
                'RowNames', obj.name, 'VariableNames', obj.name);
            
            obj.landing = array2table(zeros(ngroup, ngear), ...
                'RowNames', obj.name, 'VariableNames', obj.fleet);
            obj.discard = array2table(zeros(ngroup, ngear), ...
                'RowNames', obj.name, 'VariableNames', obj.fleet);
            
            obj.df = array2table(zeros(ngroup,ndet), ...
                'RowNames', obj.name, 'VariableNames', obj.name(nlive+1:end));
            obj.discardFate = array2table(zeros(ngear,ndet), ...
                'RowNames', obj.fleet, 'VariableNames', obj.name(nlive+1:end));
            
            svars = {'stanzaID', 'BABsplit', 'Btot', 'WmatWinf', 'RecPower'};
            
            obj.stanzadata = array2table(nan(nstan, length(svars)), ...
                'RowNames', obj.stanza, 'VariableNames', svars);
            
            warning(W);
            
            obj.pedigree = array2table(zeros(0,4), 'VariableNames', ...
                {'property', 'row', 'column', 'pedigree'});
            
        end
        
        %------------------------------------------------------------------
        function obj = set.groupdata(obj, val)
            % Check that groupdata variables meet all requirements
            
            islive = val.pp <= 1;
            
            % Prevent users from extending table
              
            if height(val) ~= obj.ngroup
                error('The groupdata table must have ngroup rows');
            end
              
            % Detritus import can only apply to detrital groups

            hasdtimp = islive & (val.dtImp > 0);
                      
            if any(hasdtimp)
                msg = sprintf('Non-zero values found for detritus import of live groups;\nvalues have been reset to 0');
                warnmessage(msg, hasdtimp, 'dtImp', val);
                val.dtImp(hasdtimp) = 0;
            end

            % No NaNs in detritus import
            
            nandtimp = islive & isnan(val.dtImp);
            if any(nandtimp)
                warnmessage('NaN found in detritus import, replacing with 0', ...
                            nandtimp, 'dtImp', val);
                val.dtImp(nandtimp) = 0;
            end
            
            % Check for negative palceholders
            
            epfields = {'b', 'pb', 'qb', 'ee', 'ge', 'bh'};
            for ifield = 1:length(epfields)
                isneg = (val.(epfields{ifield}) < 0);
                if any(isneg)
                    msg = sprintf('Negative placeholders found in %s field, replacing with NaN', epfields{ifield});
                    warnmessage(msg, isneg, epfields{ifield}, val);
                    val.(epfields{ifield})(isneg) = NaN;
                end
            end
            
            % Producers and detritus: Q/B, GE, GS should be 0
            
            qbnotzero = val.pp >= 1 & val.qb ~= 0;
            if any(qbnotzero)
                warnmessage('Non-zero value found for a producer or detrital Q/B, replacing with zero', ...
                        qbnotzero, 'qb', val);

                val.qb(qbnotzero) = 0;
            end
            
            genotzero = val.pp >= 1 & val.ge ~= 0;
            if any(genotzero)
                warnmessage('Non-zero value found for a producer or detrital GE, replacing with zero', ...
                        genotzero, 'ge', val);
                val.ge(genotzero) = 0;
            end
            
            gsnotzero = val.pp >= 1 & val.gs ~= 0;
            if any(gsnotzero)
                warnmessage('Non-zero value found for a producer or detrital GS, replacing with zero', ...
                        gsnotzero, 'gs', val);
                val.gs(gsnotzero) = 0;
            end
            
            % Detrital P/B: Should be 0 (turnover rate used in Rpath should 
            % be in detpb column)
            
            pbnotzero = val.pp == 2 & val.pb ~= 0;
            if any(pbnotzero)
                warnmessage('Non-zero value found for detrital P/B, replacing with zero (shifting to detpb column where empty)', ...
                        pbnotzero, 'pb', val);
                pbtoshift = pbnotzero & isnan(val.detpb);
                val.detpb(pbtoshift) = val.pb(pbtoshift);
                val.pb(pbnotzero) = 0;
            end
            
            % Detrital biomass doesn't effect Ecopath balance, but needs to
            % be filled in as something
            
            detbmissing = isnan(val.b) & val.pp == 2;
            if any(detbmissing)
                warnmessage('Detritus groups found with missing biomass; replacing with 0', ...
                        detbmissing, 'b', val);
                val.b(detbmissing) = 0;
            end
            
            % Area fraction should be set to something.  Assume 1 if not
            % specified.
            
            afmissing = isnan(val.areafrac);
            if any(afmissing)
                warnmessage('Missing area fraction; replacing with 1', ...
                    afmissing, 'areafrac', val);
                val.areafrac(afmissing) = 1;
            end
            
            % There are a few values in the groupdata table where you can
            % choose to set values for either one or another property
            %
            % B or BH
            % BA or BA Rate
            % Emig or Emig Rate
            %
            % Check these, make sure user doesn't set both
           
            bothb = ~isnan(val.b) & ~isnan(val.bh);
            if any(bothb)
                str = sprintf('%s, ', val.Properties.RowNames{bothb});
                str = str(1:end-2);
                error('Cannot set both B and BH (%s) for a group; choose one and set the other to NaN', str);
            end
            
            bothba = ~isnan(val.ba) & ~isnan(val.baRate);
            if any(bothba)
                str = sprintf('%s, ', val.Properties.RowNames{bothba});
                str = str(1:end-2);
                error('Cannot set both BA and BARATE (%s) for a group; choose one and set the other to NaN', str);
            end
            
            bothem = ~isnan(val.emig) & ~isnan(val.emigRate);
            if any(bothem)
                str = sprintf('%s, ', val.Properties.RowNames{bothem});
                str = str(1:end-2);
                error('Cannot set both EMIG and EMIGRATE (%s) for a group; choose one and set the other to NaN', str);
            end
            
            
%             bothb = ~isnan(val.b) & ~isnan(val.bh);
%             unsynced = val.b(bothb) ~= val.bh(bothb).*val.areafrac(bothb);
%             if any(unsynced)
%                 str = sprintf('  %s\n', val.Properties.RowNames{unsynced});
%                 error('b and bh*areafrac values are out of sync for:\n%s\nNo changes to groupdata table were made', str);
%             end
%             
%             bothe = ~isnan(val.emig) & ~isnan(val.emigRate) & ~isnan(val.b);
%             unsynced = val.emig(bothe) ~= val.emigRate(bothe).*val.b(bothe);
%             if any(unsynced)
%                 str = sprintf('  %s\n', val.Properties.RowNames{unsynced});
%                 error('emig and b*emigRate values are out of sync for:\n%s\nNo changes to groupdata table were made', str);
%             end
%             
%             bothba = ~isnan(val.ba) & ~isnan(val.baRate) & ~isnan(val.b);
%             unsynced = val.ba(bothba) ~= val.baRate(bothba).*val.b(bothba);
%             if any(unsynced)
%                 str = sprintf('  %s\n', val.Properties.RowNames{unsynced});
%                 error('ba and b*baRate values are out of sync for:\n%s\nNo changes to groupdata table were made', str);
%             end
            
            % Set groupdata after all adjustments have been made
            
            obj.groupdata = val;
            
        end
        
        %------------------------------------------------------------------
        function obj = set.pedigree(obj, val)
            if isempty(val)
                obj.pedigree = val;
            else
                idx = stanzaindices(obj);
                nonlead = cellfun(@(x) x(1:end-1), idx, 'uni', 0);
                nonlead = cat(1, nonlead{:});

                [tf, loc] = ismember({'b','pb','qb','ba','baRate','pp'}, obj.groupdata.Properties.VariableNames);
                
                % Can't apply pedigree to non-leading stanza group B, PB,
                % or QB
                
                isnonlead = strcmp(val.property, 'groupdata') & ...
                            ismember(val.column, loc(1:3)) & ...
                            ismember(val.row, nonlead);
                        
                if any(isnonlead)
                    warning('Cannot set pedigree value for non-leading stanza group B, PB, or QB; removing from table')
                end
                       
                % Can't alter stanza group BA or BARATE (for now... may
                % change as I update Rpath-related calcs), must use
                % BABsplit
                
                isstzba = strcmp(val.property, 'groupdata') & ...
                          ismember(val.column, loc(4:5)) & ...
                          ismember(val.row, cat(1, idx{:}));
               
                if any(isstzba)
                    warning('Cannot set pedigree value for stanza group BA or BARATE (use BABSPLIT in stanzadata table instead); removing from table');
                end
                
                % Can't touch stanza index
                
                issidx = strcmp(val.property, 'stanzadata') & ...
                         val.column == 1;
                     
                if any(issidx)
                    warning('Cannot set pedigree value for stanza ID; removing from table');
                end
                
                
                % Can't alter pp
                
                ispp = strcmp(val.property, 'groupdata') & ...
                              val.column == loc(6);
                if any(ispp)
                    warning('Cannot set pedigree value for pp; removing from table');
                end
                
                % Set
                
                isbad = isnonlead | isstzba | issidx | ispp;
                obj.pedigree = val(~isbad,:);
                
                % Make sure no duplicates.  If so, use last.
                
                [~, idx] = unique(obj.pedigree(:,1:3), 'last');
                obj.pedigree = obj.pedigree(idx,:);
                
                % Get rid of pedigrees assigned to a 0 or NaN
                
                vtmp = getpedigreevals(obj);
                obj.pedigree = obj.pedigree(vtmp~=0,:);      

            end
            
            
        end
        
        %------------------------------------------------------------------
        function obj = set.name(obj, val)
            validateattributes(val, {'cell'}, {'vector', 'numel', obj.ngroup});
            if ~all(cellfun(@ischar, val))
                error('Names must be cell array of string');
            end
            if ~all(cellfun(@isvarname, val))
                warning('Names must meet variable name restirctions; modfying');
                val = matlab.lang.makeValidName(val, 'delete');
            end
            
            obj.name = val;
            if ~isempty(obj.groupdata) % First set
                obj.groupdata.Properties.RowNames = val;
                obj.dc.Properties.RowNames = val;
                obj.dc.Properties.VariableNames = val;
                obj.landing.Properties.RowNames = val;
                obj.discard.Properties.RowNames = val;

                isdet = obj.groupdata.pp == 2;

                obj.df.Properties.RowNames = val;
                obj.df.Properties.VariableNames = val(isdet);
                obj.discardFate.Properties.VariableNames = val(isdet);
            end
            
        end
        
        %------------------------------------------------------------------
        function obj = set.fleet(obj, val)
            validateattributes(val, {'cell'}, {'vector', 'numel', obj.ngear});
            if ~all(cellfun(@ischar, val))
                error('Names must be cell array of string');
            end
            if ~all(cellfun(@isvarname, val))
                warning('Names must meet variable name restirctions; modfying');
                val = matlab.lang.makeValidName(val, 'delete');
            end
            
            obj.fleet = val;
            obj.landing.Properties.VariableNames = val;
            obj.discard.Properties.VariableNames = val;
            obj.discardFate.Properties.RowNames = val;
            
        end
            
        
        %------------------------------------------------------------------
        function A = checkstanza(A)
        %SETSTANZAS Fill in (or validate) B, QB, and BA values for stanzas
        %
        % A = checkstanza(A)
        %
        % This method is a wrapper around the calcstanza method.  Rather
        % than recalculating all multi-stanza group values, it checks to
        % see whether any existing non-leading-stanza values are consistent
        % with the expected values, within a small tolerance.  This allows
        % externally-calculated values to stay, even it they are off by a
        % tiny bit (can be useful when importing data from EwE6; my version
        % of the stanza calculations use a slightly different numerical
        % methods, and therefore can differ due to rounding error).
        %
        % If the values differ by more than the specified tolerance
        % (relative error > 0.005), the values are recalculated, with
        % warnings issued. 
        
            % If non-leading stanza group data was missing, fill it in (all or
            % nothing... if one stanza-group is missing, all other stanzas of that
            % group need to be recalculated)

            Tmp = A.calcstanza;

            bfill = isnan(A.groupdata.b) & ~isnan(Tmp.groupdata.b);
            bchange = ismember(A.groupdata.stanza, unique(A.groupdata.stanza(bfill)));
            A.groupdata.b(bchange) = Tmp.groupdata.b(bchange);
            qfill = isnan(A.groupdata.qb) & ~isnan(Tmp.groupdata.qb);
            qchange = ismember(A.groupdata.stanza, unique(A.groupdata.stanza(qfill)));
            A.groupdata.qb(qchange) = Tmp.groupdata.qb(qchange);

            % Check for any other changes.  Keep original data if it's within my
            % tolerance (meaning probably correct, just picking up the differences
            % between my implementation of the staza calculations vs EwE6's
            % implementation).  If it's further off, assume incorrect data, and
            % replace data for all stanzas of that group.

            berr  = (Tmp.groupdata.b  - A.groupdata.b)./A.groupdata.b;
            qerr  = (Tmp.groupdata.qb - A.groupdata.qb)./A.groupdata.qb;
            baerr = (Tmp.groupdata.ba - A.groupdata.ba)./A.groupdata.ba;

            tol = 0.005;
            bwrong = abs(berr) > tol;
            qwrong = abs(qerr) > tol;
            bawrong = abs(baerr) > tol;    
        
            warnoff = false;
            if any([bwrong; qwrong; bawrong])

                bchange = ismember(A.groupdata.stanza, unique(A.groupdata.stanza(bwrong)));

                if ~warnoff
                    msg = 'Multi-stanza data inconsistent; replacing all B, Q/B, and BA in stanza groups related to these B...';
                    tmp = [A.name(bwrong) num2cell([A.groupdata.b(bwrong) Tmp.groupdata.b(bwrong)])]';
                    str = sprintf('  %s: %.4g -> %.4g\n', tmp{:});
                    warning('%s:\n%s', msg, str);
                end
                A.groupdata.b(bchange) = Tmp.groupdata.b(bchange);

                qchange = ismember(A.groupdata.stanza, unique(A.groupdata.stanza(qwrong)));

                if ~warnoff
                    msg = '...Q/B...';
                    tmp = [A.name(qwrong) num2cell([A.groupdata.qb(qwrong) Tmp.groupdata.qb(qwrong)])]';
                    str = sprintf('  %s: %.4g -> %.4g\n', tmp{:});
                    warning('%s:\n%s', msg, str);
                end
                A.groupdata.qb(qchange) = Tmp.groupdata.qb(qchange);

                bachange = ismember(A.groupdata.stanza, unique(A.groupdata.stanza(bawrong)));
                if ~warnoff
                    msg = '... and BA';
                    tmp = [A.name(bawrong) num2cell([A.groupdata.ba(bawrong) Tmp.groupdata.ba(bawrong)])]';
                    str = sprintf('  %s: %.4g -> %.4g\n', tmp{:});
                    warning('%s:\n%s', msg, str);
                end
                A.groupdata.ba(bachange) = Tmp.groupdata.ba(bachange);

            end   
       
        end
        
        %------------------------------------------------------------------
        function idx = stanzaindices(A)
        %STANZAINDICES Extract indices of stanza groups, in order of age
        %
        % idx = stanzaindices(A)
        %
        % Input variables:
        %
        %   A:      ecopathmodel object
        %
        % Output variables:
        %
        %   idx:    nstanza x 1 cell array, each cell holds vector of group
        %           indices belonging to that stanza set, ordered from
        %           youngest to oldest (oldest  = leading)  

            grp = [(1:A.ngroup)' A.groupdata.stanza A.groupdata.ageStart];
            iss = grp(:,2) > 0;
            if any(iss)
                [sidx, grpdata] = aggregate(grp(iss,2), grp(iss,[1 3]));
                grpdata = cellfun(@(x) sortrows(x,2), grpdata, 'uni', 0);

                idx = cellfun(@(x) x(:,1), grpdata, 'uni', 0);
            else
                idx = cell(0);
            end
        end
        
        function vmid = getpedigreevals(A)
        %GETPEDIGREEVALS Extract values corresponding to pedigree entries
        %
        % vmid = getpedigreevals(A)
        %
        % Extracts the model values that correspond to each pedigree table
        % entry.
        %
        % Input variables:
        %
        %   A:      ecopathmodel object
        %
        % Output variables:
        %
        %   vmid:   nped x 1 array, where nped is the number of rows in the
        %           pedigree table
        
            nvar = height(A.pedigree);
            if nvar > 0
                [tbl, pedidx] = aggregate(A.pedigree.property, [A.pedigree.row A.pedigree.column (1:nvar)']);

                vmid = zeros(nvar,1);
                for it = 1:length(tbl)
                    tmp = table2array(A.(tbl{it}));
                    idx = sub2ind(size(tmp), pedidx{it}(:,1), pedidx{it}(:,2));
                    vmid(pedidx{it}(:,3)) = tmp(idx);
                end
            else
                vmid = [];
            end
        end


    end
end

function warnmessage(msg, mask, fld, tbl)
    tmp = [tbl.Properties.RowNames(mask) num2cell(tbl.(fld)(mask))]';
    str = sprintf('  %s (%.2f)\n', tmp{:});
    wrn = warning('off', 'backtrace');
    warning('Ecopathmodel:groupdataValidation', '%s:\n%s', msg, str);
    warning(wrn);
end
