function A = networkindices(T)
%NETWORKINDICES Calculate network indices for a food web
%
% A = networkindices(T)
%
% This function replicates the network index calculations from the
% NetIndices R Package, as documented in 
%
% Kones, J.K., Soetaert, K., van Oevelen, D. and J.Owino (2009). Are
% network indices robust indicators of food web functioning? a Monte Carlo
% approach. Ecological Modelling, 220, 370-382.  
% DOI: 10.1016/j.ecolmodel.2008.10.012
%
% Note that for this function, I expect the input matrix to group all
% imports, exports, and dissipation terms into only 3 terms, rather than
% allowing an unlimited number of out-of-system "boxes" in the food web.
% As a result, certain values that depend on the number of import/export
% nodes (such as Ltot, or A) may differ slightly from the values returned
% by the R package.
%
% Input variables:
%
%   T:      (n+3)x(n+3) flow matrix, where T(i,j) represents flow from
%           source j to sink i.  Rows/columns 1:n represent the n nodes in
%           the food web, n+1 represents the sink for flows out of the
%           system, n+2 represents the sink for respiratory flow, and n+3
%           is the source of flows into the system.
%           
%           Alternatively, T can be an (n+3)x(n+3)x(m) matrix, where the
%           third dimension represents different values for a web with the
%           same topology.  The adjacency matrix of each slice of T, i.e
%           T(:,:,m)>0, must be the same for all values of m.  Calling the
%           function with a 3D input rather than looping over many 2D
%           inputs reduces the computation time associated with certain
%           metrics that would be the same for each.
%
% Output variables:
%
%   A:      1 x 1 structure with the following fields.  Values that are
%           based only on web topology and not flux values will be scalars,
%           all others will be 1 x m arrays:
%
%           (General suite)
%
%           n:          Number of internal components
%
%           TSTp:       Total system throughput  
%
%           TSTf:       Total system throughflow
%
%           Ltot:       Number of links
%
%           Lint:       Number of internal links
%           
%           LD:         Link density
%
%           C:          Connectance
%
%           Tavg:       Average link weight
%
%           TSTfavg:    Average compartment throughflow
%
%           Cavg:       Compartmentalization
%
%           (Ascendency suite)
%
%           A:          Ascendency
%
%           DC:         Development capacity
%
%           O:          System overhead
%
%           AC:         Extent of development
%
%           (AMI suite, network uncertainty)
%
%           AMI:        Average mutual information
%
%           HR:         Statistical uncertainty
%
%           DR:         Conditional uncertainty
%
%           RU:         Realized uncertainty
%
%           Hmax:       Network uncertainty
%
%           Hc:         Constraint information
%
%           CE:         Constraint efficiency
%
%           (Finn's suite, pathway analysis)           
%
%           TSTc:       Total system cycled throughflow
%
%           TSTs:       Total system non-cycled throughflow 
%
%           FCI:        Finn's cycling index
%
%           FCIb:       Revised Finn's cycling index

% Copyright 2016 Kelly Kearney

% Check input

validateattributes(T, {'numeric'}, {'nonnegative'});

if ndims(T) > 3 || ~isequal(size(T,1),size(T,2))
    error('T must be 2D or 3D, with the first two dimensions having equal size');
end

if size(T,3) > 1
    tmp = reshape(T,[],size(T,3))>0;
    if ~(all(tmp,2) | ~any(tmp,2))
        error('T must have same adjacency matrix along the third dimension');
    end
end

% Set up output matrices

nweb = size(T,3);

ncomp = size(T,1) - 3;
export = ncomp + (1:2);
import = ncomp + 3;

flds = {'n' 'TSTp' 'TSTf' 'Ltot' 'Lint' 'LD' 'C' 'Tavg' 'TSTfavg' 'Cavg' ...
    'A' 'DC' 'O' 'AC' 'AMI' 'HR' 'DR' 'RU' 'Hmax' 'Hc' 'CE' 'Hsys' 'APL' ...
    'TSTc' 'TSTs' 'FCI' 'FCIb'};
vals = cell(length(flds),1);
[vals{:}] = deal(NaN);
issame = ismember(flds, {'n', 'Ltot', 'Lint', 'LD', 'C', 'Cavg', 'Hmax'});
[vals{~issame}] = deal(nan(1,nweb));

A = cell2struct(vals, flds);

for ii = 1:nweb

    % Setup calcs

    Tij = T(:,:,ii);
    Tint = Tij(1:ncomp,1:ncomp);

    FlowFrom  = sum(Tij,1);
    FlowTo    = sum(Tij,2);
    FlowFromC = FlowFrom(1:ncomp);
    FlowToC   = FlowTo(1:ncomp);

%     diet = bsxfun(@rdivide, Tint, FlowToC);

    %------------------------
    % General network indices
    %------------------------

    if ii == 1
        A.n = ncomp;                    % Number of internal compartments
    end
    A.TSTp(ii) = sum(Tij(:));           % Total system throughput

    ratecomp = FlowToC' - FlowFromC; 

    % Total system throughflow

    A.TSTf(ii) = sum(Tint(:)) + sum(FlowFrom(import)) - sum(ratecomp(ratecomp < 0));

    if ii == 1
        A.Ltot = nnz(Tij);              % Number of links,
        A.Lint = nnz(Tint);             % Number of internal links
        A.LD = A.Ltot./A.n;             % Link density
        A.C = A.Lint./(A.n*(A.n-1));    % Connectance
    end

    A.Tavg(ii) = A.TSTp(ii)./A.Ltot;    % Average link weight
    A.TSTfavg(ii) = A.TSTf(ii)./A.n;    % Average compartment throughflow

    % Compartmentalization

    if ii == 1
        adj = Tint > 0;
        c = zeros(size(adj));
        for ir = 1:A.n
            for jj = 1:A.n
                pp1 = adj(:,ir) | adj(ir,:)'; % groups #1 interacts with
                pp2 = adj(:,jj) | adj(jj,:)'; % groups #2 interacts with

                c(ir,jj) = sum(pp1 & pp2)./sum(pp1 | pp2);
            end
        end

        A.Cavg = (sum(c(:)) - A.n)./(A.n*(A.n-1));
    end

    %------------------------
    % Growth and development
    % Indices
    %------------------------

    tmp = bsxfun(@times, FlowFrom, FlowTo);

    asc = Tij.*log2(Tij.*A.TSTp(ii)./tmp);
    cap = Tij.*log2(Tij./A.TSTp(ii));
    ispos = Tij > 0;

    A.A(ii) = sum(asc(ispos));      % Ascendency
    A.DC(ii) = -sum(cap(ispos));    % Development capacity
    A.O(ii) = A.DC(ii) - A.A(ii);   % Overhead
    A.AC(ii) = A.A(ii)./A.DC(ii);   % Extent of development

    %------------------------
    % Network uncertainty 
    % indices
    %------------------------

    A.AMI(ii) = A.A(ii)/A.TSTp(ii); % Average mutual information

    Q = FlowFrom./A.TSTp(ii);
    Q = Q(Q > 0);

    A.HR(ii) = -sum(Q.*log2(Q));    % Statistical uncertainty
    A.DR(ii) = A.HR(ii) - A.AMI(ii);% Conditional uncertainty index
    A.RU(ii) = A.AMI(ii)/A.HR(ii);  % Realized uncertainty index
    A.Hmax = ncomp * log2(ncomp+1); % Network uncertainty

    blsum = bsxfun(@rdivide, Tij, FlowFrom) .* log2(bsxfun(@rdivide, Tij, FlowFrom));
    blsum = sum(blsum(Tij>0));

    A.Hc(ii) = blsum + A.Hmax;      % Constraint information
    A.CE(ii) = A.Hc(ii)./A.Hmax;    % Constraint efficiency
    A.Hsys(ii) = A.Hmax - A.Hc(ii); % Network efficiency

    %------------------------
    % Pathway analysis
    % indices
    %------------------------

    % Average path length

    A.APL(ii) = A.TSTf(ii)./ (sum(FlowTo(export)) + sum(ratecomp(ratecomp>0)));

    comptf = max(FlowFromC', FlowToC);

    Qij = bsxfun(@rdivide, Tint, comptf);
    I = eye(size(Tint));
    IQ = I - Qij;
    M = inv(IQ);
    diaM = diag(M);

    A.TSTc(ii) = sum((1-1./diaM).*FlowFromC');  % Cycled throughflow
    A.TSTs(ii) = A.TSTf(ii) - A.TSTc(ii);       % Non-cycled throughflow
    A.FCI(ii) = A.TSTc(ii)/A.TSTf(ii);          % Finn's Cycling Index
    A.FCIb(ii) = A.TSTc(ii)/A.TSTp(ii);         % Finn's Cycling Index, revised
end

