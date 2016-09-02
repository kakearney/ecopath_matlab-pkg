function G = graph(EM, varargin)
%GRAPH Convert ecopathmodel object to a digraph object
%
% G = graph(EM)
% G = graph(EM, p1, v1, ...)
%
% Input variables:
%
%   EM:     ecopathmodel object
%
% Optional input parameters (passed as parameter/value pairs, defaults in
% []): 
%
%   oos:    logical scalar, true to include out-of-system nodes [true]
%
%   det:    logical scalar, true to include flow-to-detritus edges [true]
%
% Output variables:
%
%   G:      digraph object.  Nodes represent the functional groups from the
%           Ecopath model, plus 4 out-of-system nodes: Import/Export
%           (source/target for import/export), Respiration (target for
%           respiration losses from living groups), Landings (target for
%           fisheries landings removed from the system), and Primary
%           production (source for primary production).  The Nodes table
%           includes the following fields:
%
%           Name:   name of group
%           B:      biomass (same mass units as input model, typically t ww
%                   km^-2) 
%           type:   0 = consumer
%                   1 = producer
%                   2 = detritus
%                   3 = fleet
%                   4 = out-of-system
%           TL:     trophic level (note: same as ecopath-calculated value
%                   for living/detrital groups, but also adds fishing
%                   fleets to the calculation, treating them as a
%                   predator of both their landed and discarded catches).
%
%           Edges represent the fluxes between nodes, including predation,
%           primary production, respiration, flows to detritus,
%           import/export, fisheries catches, fisheries landings, and
%           fisheries discards.  The Edges table includes the following
%           fields:
%
%           EndNodes:   source and target nodes for each flux
%           Weight:     mass flux per unit time (same units as input model,
%                       typically t ww km^-2 yr^-1) 
%

% Copyright 2016 Kelly Kearney

% Main nodes: groups, fleets, and out-of-system sources and sinks 

p = inputParser;
p.addParameter('oos', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('det', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.parse(varargin{:});
Opt = p.Results;

Ep = EM.ecopath;

outofsys = {'Import/export', 'Respiration', 'Landings', 'Primary production'}';
noos = length(outofsys);

nodenames = [EM.name; EM.fleet; outofsys];
nnode = length(nodenames);

G = digraph(Ep.flow, nodenames);

% Add a few properties from Ecopath calculations
% Type: 0 = consumer, 1 = producer, 2 = detritus, 3 = fleet, 
%       4 = out-of-system

bfleet = sum(table2array(EM.landing), 1)'; % Don't really have a biomass, so counting "kept" biomass

G.Nodes.B = [Ep.b; bfleet; zeros(noos,1)];
G.Nodes.type = [ceil(EM.groupdata.pp); ones(EM.ngear,1)*3; ones(noos,1)*4]; 

% Redo trophic level calculations, including groups and gear

adj = zeros(EM.ngroup+EM.ngear);
adj(1:EM.ngroup,1:EM.nlive) = Ep.q0(1:EM.ngroup,1:EM.nlive);
adj(1:EM.ngroup,(1:EM.ngear)+EM.ngroup) = table2array(EM.landing)+table2array(EM.discard);
dc = bsxfun(@rdivide, adj, sum(adj,1));    
dc(isnan(dc)) = 0;

pp = [EM.groupdata.pp; zeros(EM.ngear,1)];

[pp, isrt] = sort(pp);
dc = dc(isrt,isrt);

tl = trophiclevel(dc, pp, EM.nlive+EM.ngear, EM.ngroup+EM.ngear);
[~,loc] = ismember(1:(EM.ngroup+EM.ngear), isrt);
G.Nodes.TL = [tl(loc); nan(noos,1)];

if ~Opt.oos
    G = rmnode(G, outofsys);
end

if ~Opt.det
    detnode = G.Nodes.Name(G.Nodes.type == 2);
    isdet = ismember(G.Edges.EndNodes(:,2), detnode);
    G = rmedge(G, find(isdet));
end
