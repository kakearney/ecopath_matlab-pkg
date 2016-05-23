classdef ecopathmodel
%ECOPATHMODEL A Matlab-based version of the Ecopath food web model
%   
%
    
    properties
        ngroup;
        %NGROUP Number of groups (living and detrital) in model
        nlive;
        ngear;
        name;
        fleet;
        groupdata;
        dc;
        landing;
        discard;
        df;
        discardFate;
        stanza;
        stanzadata;
    end
    
    methods
        
        % Constructor
        function obj = ecopathmodel(ngroup, nlive, ngear, varargin)
            
            grpdefault = strtrim(cellstr(num2str((1:ngroup)', 'group%d')));
            fltdefault = strtrim(cellstr(num2str((1:ngear)', 'fleet%d')));
            stzdefault = {};

            p = inputParser;
            addParameter(p, 'groups', grpdefault, @(x) iscell(x) && length(x)==ngroup && all(cellfun(@ischar, x)));
            addParameter(p, 'fleets', fltdefault, @(x) iscell(x) && length(x)==ngear && all(cellfun(@ischar, x)));
            addParameter(p, 'stanzas',stzdefault); 
            parse(p, varargin{:});
            
            obj.ngroup = ngroup;
            obj.nlive = nlive;
            obj.ngear = ngear;
            
            gvars = {'b', 'pb', 'qb', 'ee', 'ge', 'gs', 'dtImp', 'bh', ...
                    'pp', 'areafrac', 'ba', 'baRate', 'immig', 'emig', ...
                    'emigRate', 'stanza', 'ageStart', 'vbK', 'detpb'};
            nv = length(gvars);
            
            ndet = ngroup - nlive;
            
            obj.name = p.Results.groups;
            obj.fleet = p.Results.fleets;
            obj.stanza = p.Results.stanzas;
            
            nstan = length(obj.stanza);
            
            obj.groupdata = array2table(nan(ngroup,nv), ...
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
            
            svars = {'stanzaID', 'BABsplit'};
            
            obj.stanzadata = array2table(nan(nstan, length(svars)), ...
                'RowNames', obj.stanza, 'VariableNames', svars);
            
        end
        
      
        
    end
    
end

