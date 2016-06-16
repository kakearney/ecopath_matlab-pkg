function A = calcstanza(A, varargin)
% CALCSTANZA Calculate B and Q/B values for multi-stanza Ecopath groups
%
% A = calcstanza(EM)
% A = calcstanza(EM, p1, v1)
%
% This function calculates the B and QB values associated with multistanza
% ecopath groups.
%
% A few caveats...
%
% The values calculated by this script do not exactly replicate the values
% one will get using the EwE 6 "Edit multi-stanza" menu.  This is the
% result of both differing precisions between EwE6 and this code (I use
% double, EwE6 uses a mix of integer and single) and a slightly different
% method of dealing with the upper tail end of the biomass/consumption rate
% age curves (I extend to 99.99% relative Winf; EwE6 extends to 90% and
% then adds an extra accumulator calculation to deal with the rest). In
% tests, the percent error between my results and the EwE ones tend to be
% <0.01% for biomass values, and 0.005-0.4% for consumption rate values.
% While the differences are very small, they may be important if one plans
% to generate ensembles using this code and then move back to the EwE
% software. If you do this, I suggest you double-check the balance of the
% model(s) there. 
%
% Input variables:
%
%   EM:     ecopathmodel object 
%
% Optional input variables (passed as parameter/value pairs)
%
%   plot:   logical scalar, true to plot growth curve details for each
%           stanza group, primarily for debugging purposes.  [false] 
%
%   da:     discretizatin interval (months) for age calculations [1]
%
% Output variables:
%
%   A:      ecopathmodel object, identical to EM except that
%           non-leading multi-stanza group B and QB values have been
%           recalculated/filled in.  

% Copyright 2016 Kelly Kearney

% Parse input

p = inputParser;
p.addParameter('plot', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('da', 1, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.parse(varargin{:});

Opt = p.Results;

ns = size(A.stanzadata,1);

% Set up plotting (if true)

if Opt.plot
    
    h.fig = figure;
    
    vfrac = 0.9;
    hght1 = vfrac./ns;
    hght2 = hght1*0.8;
    
    yax = 0.05:hght1:0.9;
    yax = yax(end:-1:1);
    
    if verLessThan('matlab', 'R2014b')
        h.ax = zeros(ns,1);
    else
        h.ax = gobjects(ns,1);
    end
    for ii = 1:ns
        h.ax(ii) = axes('position', [0.05 yax(ii) 0.9 hght2]);
    end
end

% Make sure all multi-stanza groups have B and BA (not BH or BARATE) set. 
% TODO: Right now BAsplit override BA... still need to figure this out.

calcb = A.groupdata.stanza > 0 & isnan(A.groupdata.b) & ~isnan(A.groupdata.bh);
btmp = A.groupdata.areafrac(calcb) .*  A.groupdata.bh(calcb);
A.groupdata.bh(calcb) = NaN;
A.groupdata.b(calcb) = btmp;

% Loop over multi-stanza groups

for is = 1:ns
    
    idx = find(A.groupdata.stanza == is);
    
    [a, isrt] = sort(A.groupdata.ageStart(idx));
    idx = idx(isrt);
   
    % Parameters
    
    k = A.groupdata.vbK(idx(1));               % Curvature parameter
    if any(isnan(A.stanzadata.BABsplit))
        bab = A.groupdata.baRate(idx);
        isn = isnan(bab);
        bab(isn) = A.groupdata.ba(idx(isn))./A.groupdata.b(idx(isn));
    else
        bab = A.stanzadata.BABsplit(is);
    end
    z = A.groupdata.pb(idx);
    if isnan(A.groupdata.b(idx(end)))
        blead = A.stanzadata.Btot(is);
        bflag = 'total';
    else
        blead = A.groupdata.b(idx(end));
        bflag = 'lead';
    end
    qblead = A.groupdata.qb(idx(end));
    
    % Fill in non-leading group values
    
    [Out,D] = editstanzacalcs(a, k, bab, blead, qblead, z, Opt.da, bflag);
    
    A.groupdata.b(idx)  = Out.b;
    A.groupdata.qb(idx) = Out.qb;
    A.groupdata.baRate(idx) = NaN; % Remove BARATE if present, going to use BA instead
    A.groupdata.ba(idx) = Out.ba;
    
    % Plot to check, mimicking the one in EwE6
    
    if Opt.plot
        
        yfac = max(D.y);
        ynorm = bsxfun(@rdivide, D.y, yfac);
        line(D.x, ynorm, 'parent', h.ax(is));
        line([a a]', (ones(size(a))*[0 1])', 'color', 'k', 'parent', h.ax(is));
        text(max(D.x)*0.98, 0.5, A.stanza{is}, ...
            'parent', h.ax(is), 'vert', 'top', 'horiz', 'right', ...
            'fontsize', 8);
        set(h.ax(is), 'xlim', [0 max(D.x)]);     
    end
    
end

% Add labels to plots

if Opt.plot
    xlabel(h.ax(end), 'Age (months)');
    lbl = {...
            'w = von Bertalanffy body weight'
            'l = survivorship'
            'w*l'
            '$\sum_{0}^{a}{Z_a} - a\frac{BA}{B}$'};
    hleg = legendflex(h.ax(1), lbl, 'ref', h.ax(1), 'nrow', 1, 'anchor', {'n','s'}, ...
        'xscale', 0.5, 'buffer', [0 5], 'interpreter', 'latex');
%     set(findall(hleg, 'type', 'text'), 'interpreter', 'Latex');
end
    






