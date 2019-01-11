function S = networkindices(EM, varargin)
%NETWORKINDICES Calculate ecological network indices
%
% S = networkindices(EM)
% S = networkindices(EM, p1, v1, ...)
%
% This function calculates a variety of network indices for an Ecopath
% model.  It was inspired by the following study:
%
% Guesnet V, Lassalle G, Chaalali A, Kearney K, Saint-B�at B, Karimi B,
% Grami B, Tecchio S, Niquil N, Lobry J (2015) Incorporating food-web
% parameter uncertainty into Ecopath-derived ecological network indicators.
% Ecol Modell 313:29?40   
% DOI: 10.1016/j.ecolmodel.2015.05.036
%
% The specific network indices calculated herein have been streamlined
% compared to the code offered alongside the above publication. The indices
% calculated are primarily based on the NetIndices R package; these are
% similar to (but in some cases not identical to) the formulas used by
% Statistics and Network Analysis tools in EwE6.  I have also included a
% handful of Ecopath-specific indices that are not part of the NetIndices
% package.    
%
% Input variables:
%
%   EM:         ecopathmodel object
%
% Optional input variables:
%
%   fleet:      How are fishing fleets treated?
%               'in':   Considered nodes of the system, catches are
%                       within-system fluxes, landings are considered
%                       exports from the fleet node and discards to
%                       detritus are within-system fluxes  
%               'out':  Fleets are external to the system, the landing
%                       portion of catches is an export and the discard
%                       portion is routed directly from fish group to
%                       detritus
%
%   ensemble:   nped x nset array of values.  Each column represents one
%               set of values to be substituted into EM using the
%               subpedigreevalues method.  If included, statistics will be
%               calculated for all ensemble members rather than the central
%               model defined by the ecopathmodel object.
%
% Output variables:
%
%   S:          n x 1 structure including all fields in the networkindices
%               function (the non-ecopathmodel method) output, as well as
%               the following additional ones: 
%
%               Qsum:       Sum of consumption
%
%               Psum:       Sum of production
%
%               catchTL:    mean trophic level of catch
%
%               GE:         gross efficiency (catches/net primary
%                           production) 
%
%               T:          ng x ng x nset array, Flow from compartment j
%                           to i, where j represents the columns of the
%                           flow matrix and i the rows (note that this
%                           row/column convention is the opposite of most
%                           of my Ecopath-related functions, but is
%                           consistent with the NetIndices package on which
%                           the main calculations are based).  The
%                           rows/columns correspond to in-system nodes
%                           (functional groups + gear for fleet='in',
%                           functional groups only for fleet='out'),
%                           export sink, dissipation sink, and input
%                           source, respectively.  

% Copyright 2016 Kelly Kearney

p = inputParser;
p.addParameter('fleet', 'in');
p.addParameter('ensemble', zeros(height(EM.pedigree),0), @(x) validateattributes(x, {'numeric'}, {'2d', 'nrows', height(EM.pedigree)})); 
p.parse(varargin{:});

Opt = p.Results;
Opt.fleet = validatestring(Opt.fleet, {'in', 'out'});

if isempty(Opt.ensemble)
    Ep = EM.ecopath;
    B = struct;
else
    [~, Ep] = EM.ecopath('ensemble', Opt.ensemble);
    B = EM.subpedigreevalues(Opt.ensemble);
end

nens = length(Ep);

% Set up tranport matrix for each ensemble member

for ii = nens:-1:1

    flds = {'landing', 'discard', 'discardFate'};
    for ifld = 1:length(flds)
        if isfield(B, flds{ifld})
            In.(flds{ifld}) = B.(flds{ifld})(:,:,ii);
        else
            In.(flds{ifld}) = table2array(EM.(flds{ifld}));
        end
    end

    switch Opt.fleet
        case 'in'
            if ii == nens
                exidx = [Ep(ii).Idx.out Ep(ii).Idx.lan];
                inidx = [Ep(ii).Idx.gpp];
                uuidx = Ep(ii).Idx.res;

                idx = setdiff(1:size(Ep(ii).flow,1), [exidx inidx uuidx]);
                n = length(idx);
                T = zeros(n+3,n+3,nens);
              
            end

            T(1:n,1:n,ii) = Ep(ii).flow(idx,idx);

            T(1:n,n+1,ii) = sum(Ep(ii).flow(idx,exidx), 2);
            T(1:n,n+2,ii) = Ep(ii).flow(idx,uuidx);
            T(n+3,1:n,ii) = Ep(ii).flow(inidx,idx);

        case 'out'

            if ii == nens
                exidx = [Ep(ii).Idx.out Ep(ii).Idx.lan];
                inidx = [Ep(ii).Idx.gpp];
                uuidx = Ep(ii).Idx.res;

                idx = setdiff(1:size(Ep(ii).flow,1), [exidx inidx uuidx Ep(ii).Idx.gear]);
                n = length(idx);
                T = zeros(n+3,n+3,nens);
            end
            
            % First, reroute flows so bypass gears

            catches = [sum(In.landing,2) In.discard];
            fate1 = [zeros(1,EM.ngroup-EM.nlive); In.discardFate]; 
            fate1 = [fate1 1-sum(fate1,2)];

            ttmp = Ep(ii).flow;
            ttmp(Ep(ii).Idx.gear,:) = 0;
            ttmp(:,Ep(ii).Idx.gear) = 0;
            ttmp([Ep(ii).Idx.liv Ep(ii).Idx.det], [Ep(ii).Idx.det, Ep(ii).Idx.lan]) = ...
                ttmp([Ep(ii).Idx.liv Ep(ii).Idx.det], [Ep(ii).Idx.det, Ep(ii).Idx.lan]) + catches*fate1;

            % Convert to T

            
            T(1:n,1:n,ii) = ttmp(idx,idx);

            T(1:n,n+1,ii) = sum(ttmp(idx,exidx), 2);
            T(1:n,n+2,ii) = ttmp(idx,uuidx);
            T(n+3,1:n,ii) = ttmp(inidx,idx);
    end
end

T = permute(T, [2 1 3]); % T_ij, i = sink, j = source

if any(T(:) < 0)
    error('Negative fluxes found in flow matrix T');
end

% Calculate main indices
    
S = networkindices(T);
    
% Add a few Ecopath-specific indices

S.Qsum = nan(1,nens);
S.Psum = nan(1,nens);
S.catchTL = nan(1,nens);
S.GE = nan(1,nens);

S.T = T;

for ii = 1:nens
    S.Qsum(ii) = sum(Ep(ii).q0sum(1:EM.nlive)); % Sum of consumption
    S.Psum(ii) = sum(Ep(ii).pb .* Ep(ii).b); % Sum of production
    
    fmort = Ep(ii).fishMortRate .* Ep(ii).b;
    S.catchTL(ii) = sum(fmort./sum(fmort) .* Ep(ii).trophic);
    
    netpp = sum(Ep(ii).pb .* Ep(ii).b .* (EM.groupdata.pp == 1));
    S.GE(ii) = sum(fmort)./netpp;
end





