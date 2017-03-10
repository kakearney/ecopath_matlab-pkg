function [x, Out] = createensemble(A, nsample, varargin)
%CREATEENSEMBLE Build an ensemble of Ecopath model parameters
%
% [x, Out] = createensemble(EM, nset)
% [x, Out] = createensemble(EM, nset, p1, v1, ...)
%
% This method generates an ensemble of Ecopath model parameters.  It was
% designed as a way to explore the parameter space of an Ecopath model,
% given the often-large uncertainty associated with its input parameters.
%
% This concept is similar to the first step of the Ecoranger routine
% available in EwE, but this function makes no attempt to define any set of
% parameters as the "best."
%
% Uncertainty ranges are applied to all parameters listed in the pedigree
% table of the input ecopathmodel object.
%
% Input parameters:
%
%   EM:         an ecopathmodel object.
%
%   nset:       number of ensemble members (i.e. sets of parameters) to
%               generate.
%
% Optional input parameters (passed as parameter/value pairs, defaults in
% []): 
%
%   collect:    'all' or 'balanced', specifies whether to return all
%               parameter sets, or only those that meet the Ecopath mass
%               balance criteria.  In the latter case, parameter sets will
%               be generated until nset balanced sets are found; this may
%               take a while, depending on the particulars of the input
%               ecosystem. ['balanced']     
%
%   pdfname:    Specifies the probability distribution from which
%               parameters will be chosen.  Note that this distribution
%               applied to the values fed into the subpedigreevalues
%               method; values returned after adjustments (diet
%               normalization, stable age curve calculations, etc.) may no
%               longer conform to these distributions. ['uniform']
%
%               'uniform':      Uniform distribution between x-ped*x and
%                               x+ped*x, where x is the point value from
%                               the original Ecopath model and ped is the
%                               pedigree value.
%
%               'lognormal':    Lognormal distribution with mean x and
%                               variance (ped*x/2)^2. 
%
%               'triangular':   triangular distribution with a peak at x
%                               and decreasing to 0 at x-ped*x and x+ped*x.
%
%   sample:     Specifies the technique used to choose samples from the
%               parameter distribution functions. ['mcs']
%
%               'mcs':      Monte Carlo simulation
%
%               'lhs':      Latin hypercube simulation, with a maximin
%                           criterion to maximize point distance between
%                           ensemble members
%
%   lhsiter:    For latin hypercube sampling only, number of iterations to
%               perform per sample block in an attempt to minimize distance
%               between points (i.e. 'iterations' parameter to pass to
%               lhsdeign.m) [20]
%
%   testsize:   For collect = 'balanced' only, number of ensemble parameter
%               sets to test for balance at one time.  Testing larger sets
%               is more efficient, since each time the ecopath method is
%               called, some setup calculations need to be repeated.
%               Increasing the default is particularly recommended when
%               nset is small but the input model is precariously balanced
%               (i.e. the vast majority of randomly-chosen parameter sets
%               will fail the mass balance criteria). Decreasing the
%               default may be useful if you want to increase the frequency
%               of updates to the progress-tracker printed to the screen.
%               [nsample]   
%
%   maxiter:    For collect = 'balanced' only, maximum number of sample
%               sets to test for balance (at maximum, testsize*maxiter sets
%               of parameters tested).  If provided, the requested number
%               of sets may not be found, and the output x will only have
%               as many rows as found when the limit is reached.
%
% Output variables:
%
%   x:          nset x nped array of parameter values.  Each row represents
%               one set of parameters that can be substituted into the
%               input ecopathmodel object (see subpedigreevalues method),
%               columns correspond to rows in the pedigree table of the
%               input ecopathmodel object.
%
%   Out:        1 x 1 structure with the following fields:
%
%               vmid:   nset x 1 array, point estimate (from center model)
%                       for each parameter
%
%               mu:     nset x 1 array, mu values associated with the
%                       lognormal distribution for each parameter
%
%               sig:    nset x 1 array, sigma values associated with the
%                       lognormal distribution of each parameter 
%
%               lo:     nset x 1 array, the lower-end cutoff for the
%                       uniform  and triangular distributions of each
%                       parameter
%
%               hi:     nset x 1 array, the upper-end cutoff for the
%                       uniform and triangular distributions of each
%                       parameter 
%
%               nall:   scalar, number of total parameter sets generated in
%                       order to create the required number of parameter
%                       sets.  For collect = 'all', this will be the same
%                       as nset; for collect = 'balanced', it may be much
%                       higher.

% Copyright 2016 Kelly Kearney

% Parse input

p = inputParser;
p.addParameter('pdfname', 'uniform');
p.addParameter('collect', 'balanced');
p.addParameter('maxiter', Inf, @(x) validateattributes(x, {'numeric'}, {'scalar','integer'}));
p.addParameter('sample', 'mcs');
p.addParameter('lhsiter', 20, @(x) validateattributes(x, {'numeric'}, {'scalar','integer'}));
p.addParameter('testsize', nsample, @(x) validateattributes(x, {'numeric'}, {'scalar','integer'}));

p.parse(varargin{:});

Opt = p.Results;

Opt.pdfname = validatestring(Opt.pdfname, {'uniform', 'lognormal', 'triangular'}, 'createensemble', 'pdfname');
Opt.collect = validatestring(Opt.collect, {'all', 'balanced'}, 'createensemble', 'collect');
Opt.sample = validatestring(Opt.sample, {'mcs', 'lhs'}, 'createensemble', 'sample');

%----------------
% Gather vars
%----------------

% Get values corresponding to each pedigree entry

vmid = getpedigreevals(A);
nvar = length(vmid);

vped = A.pedigree.pedigree;
vciv = vmid.*vped;

% Convert to mu/sigma values (for lognormal) ...

vvar = (vciv./2).^2; % variance

mu = log((vmid.^2)./sqrt(vvar+vmid.^2));
sigma = sqrt(log(vvar./(vmid.^2)+1));

% ... or upper/lower (for uniform, triangular)

lo = vmid - vciv;
hi = vmid + vciv;


%----------------
% Generate sets
%----------------

switch Opt.collect
    
    case 'all'
        
        % Generate samples
        
        p = samplepoints(Opt.sample, nsample, nvar, Opt.lhsiter);
        x = transformpoints(Opt.pdfname, p, lo, hi, mu, sigma, vmid);

        nall = nsample;
        
    case 'balanced'
        
        nbal = 0;
        nall = 0;
        x = zeros(nsample, nvar);
        
        [~,gscol] = ismember('gs', A.groupdata.Properties.VariableNames);
        isgs = strcmp(A.pedigree.property, 'groupdata') & A.pedigree.column == gscol;
        if any(isgs)
            gstmp = repmat(A.groupdata.gs, 1, nsample);
        end
            
        
        cpb = ConsoleProgressBar();
        cpb.setMaximum(nsample);
        fprintf('Generating ensembles...\n');
        cpb.start();

        count = 1;
        while nbal < nsample && count <= Opt.maxiter
            
            % Generate samples

            ptmp = samplepoints(Opt.sample, Opt.testsize, nvar, Opt.lhsiter);
            xtmp = transformpoints(Opt.pdfname, ptmp, lo, hi, mu, sigma, vmid);

            % Keep only balanced sets
            
            [~, Ep] = A.ecopath('ensemble', xtmp', 'skipextra', true);
            
            if any(isgs)
                gstmp(A.pedigree.row(isgs),:) = xtmp(:,isgs)';
                isbal = all([Ep.ee] <= 1 & [Ep.ee] >= 0 & ~isnan([Ep.ee]) & ...
                            [Ep.ge] >= 0 & [Ep.ge] <= 1 & ...
                            gstmp + [Ep.ge] <= 1, 1);
                
            else
                isbal = all([Ep.ee] <= 1 & [Ep.ee] >= 0 & ~isnan([Ep.ee]) & ...
                            [Ep.ge] >= 0 & [Ep.ge] <= 1 & ...
                            bsxfun(@plus, A.groupdata.gs, [Ep.ge]) <= 1, 1);
            end
            
            xbal = xtmp(isbal,:);
            nnew = size(xbal,1);
            
            if nnew + nbal <= nsample
                x((1:nnew)+nbal,:) = xbal;
                nall = nall + Opt.testsize;
            else
                x((nbal+1):end,:) = xbal(1:(nsample-nbal),:);
                nall = nall + (nsample-nbal);
            end
            
            nbal = nbal + nnew;
            
            cpb.setValue(min(nbal,nsample));
            cpb.setText(sprintf('Attempted: %d', nall));
            
            count = count + 1;
            
        end
        cpb.stop();
        fprintf('\n');    
        
        if nbal < nsample
            x = x(:,1:nbal);
        end
        
end

if nargin > 1
    Out.vmid = vmid;
    Out.mu = mu;
    Out.sig = sigma;
    Out.lo = lo;
    Out.hi = hi;
    Out.nall = nall;
end

%--------------------
% Subfunctions
%--------------------

% Generate prescribed number of uniformly-distributed points

function p = samplepoints(method, nsample, nvar, niter)

switch method
    case 'mcs'
       p = rand(nsample, nvar);
    case 'lhs'
       p = lhsdesign(nsample, nvar, ...
           'criterion', 'maximin', ...
           'iterations', niter);
end


% Tranform uniformly-distributed points to the prescribed PDF

function x = transformpoints(pdfname, p, lo, hi, mu, sigma, vmid)

ismissing = isnan(vmid);

switch pdfname
            
    case 'uniform'
        x = zeros(size(p));
        for iv = 1:size(p,2)
            x(:,iv) = unifinv(p(:,iv), lo(iv), hi(iv));
        end

    case 'lognormal'
        x = zeros(size(p));
        for iv = 1:size(p,2)
            x(:,iv) = logninv(p(:,iv), mu(iv), sigma(iv));
        end

    case 'triangular'
        x = nan(size(p));
        for iv = 1:size(p,2)
            if ~isnan(vmid(iv))
                x(:,iv) = icdf('Triangular', p(:,iv), lo(iv), vmid(iv), hi(iv));
            end
        end
end
x(:,ismissing) = NaN;


