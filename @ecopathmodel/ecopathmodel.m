classdef ecopathmodel
%ECOPATHMODEL A Matlab-based version of the Ecopath food web model
%   
%
    
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
        %               last.    
        %
        %   fleets:     ngear x 1 cell array of strings, names
        %               corresponding to each fishing gear. 
        %
        %   stanzas:    nstanza x 1 cell array of strings, where nstanza
        %               refers to the number of multi-stanza sets of
        %               groups.  Functional groups that are part of a
        %               multi-stanza set have certain constraints set on
        %               their B, P/B, and Q/B values.  
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
            
            svars = {'stanzaID', 'BABsplit', 'Btot'};
            
            obj.stanzadata = array2table(nan(nstan, length(svars)), ...
                'RowNames', obj.stanza, 'VariableNames', svars);
            
            warning(W);
            
            obj.pedigree = array2table(zeros(0,4), 'VariableNames', ...
                {'table', 'row', 'column', 'value'});
            
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
            
            % Check the either/or properties, make sure they aren't set out
            % of sync
            
            bothb = ~isnan(val.b) & ~isnan(val.bh);
            unsynced = val.b(bothb) ~= val.bh(bothb).*val.areafrac(bothb);
            if any(unsynced)
                str = sprintf('  %s\n', val.Properties.RowNames{unsynced});
                error('b and bh*areafrac values are out of sync for:\n%s\nNo changes to groupdata table were made', str);
            end
            
            bothe = ~isnan(val.emig) & ~isnan(val.emigRate) & ~isnan(val.b);
            unsynced = val.emig(bothe) ~= val.emigRate(bothe).*val.b(bothe);
            if any(unsynced)
                str = sprintf('  %s\n', val.Properties.RowNames{unsynced});
                error('emig and b*emigRate values are out of sync for:\n%s\nNo changes to groupdata table were made', str);
            end
            
            bothba = ~isnan(val.ba) & ~isnan(val.baRate) & ~isnan(val.b);
            unsynced = val.ba(bothba) ~= val.baRate(bothba).*val.b(bothba);
            if any(unsynced)
                str = sprintf('  %s\n', val.Properties.RowNames{unsynced});
                error('ba and b*baRate values are out of sync for:\n%s\nNo changes to groupdata table were made', str);
            end
            
            % Set groupdata after all adjustments have been made
            
            obj.groupdata = val;
            
        end
        
        %------------------------------------------------------------------
        function A = setstanzas(A)
        %SETSTANZAS Fill in (or validate) B, QB, and BA values for stanzas
        %
        % A = setstanzas(A)
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
        % (0.005), the values are recalculated, with warnings issued.
        
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
                    tmp = [A.name(bawrong) num2cell([A.groupdata.ba(bawrong) Tm.groupdatap.ba(bawrong)])]';
                    str = sprintf('  %s: %.4g -> %.4g\n', tmp{:});
                    warning('%s:\n%s', msg, str);
                end
                A.groupdata.ba(bachange) = Tmp.groupdata.ba(bachange);

            end   
       
        end
    end
end

function warnmessage(msg, mask, fld, tbl)
    tmp = [tbl.Properties.RowNames(mask) num2cell(tbl.(fld)(mask))]';
    str = sprintf('  %s (%.2f)\n', tmp{:});
    warning('Ecopathmodel:groupdataValidation', '%s:\n%s', msg, str);
end
