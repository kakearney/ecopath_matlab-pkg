function displaybasic(In, Ep, cflag)
%DISPLAYBASIC Prints ecopath results to screen
%
% displayecopath(EM)
% displayecopath(EM, Ep)
% displayecopath(EM, Ep, cflag)
%
% This function displays main parameters of an ecopathmodel in the
% command window, using a table structure similar to the way Ecopath
% results are presented in the original EwE "Basic Estimates" gui.
%
% Input variables:
%
%   EM:     ecopathmodel object
% 
%   Ep:     ecopath results structure (see ecopath method for details).
%           If not included, values set to NaN in the input structure will
%           be displayed as empty. If included, the calculated values are
%           shown in blue; values indicative of an unbalanced model are
%           shown in red.
%
%   cflag:  flag to print using color (via cprintf function).  The
%           cprintf function uses tags specifically formatted for the
%           Matlab command window; if printing to a different environment
%           (e.g., publishing a document, running in nodesktop mode), you
%           may want to turn this off.  Default is true.

% Copyright 2016 Kelly Kearney

% Parse input

narginchk(1,3);

epflag = nargin > 1 && ~isempty(Ep);
if nargin < 3
    cflag = true;
end


% Display basic input/estimates (TL HA BH B PB QB EE GE GS DI)

blank = nan(In.ngroup,1);

basicin = [blank ...
           In.groupdata.areafrac ...
           blank ...
           In.groupdata.b ...
           In.groupdata.pb ...
           In.groupdata.qb ...
           In.groupdata.ee ...
           In.groupdata.ge ...
           In.groupdata.gs ...
           In.groupdata.dtImp];
       
if epflag
    basicout = [Ep.trophic ...
               Ep.areafrac ...
               Ep.bh ...
               Ep.b ...
               Ep.pb ...
               Ep.qb ...
               Ep.ee ...
               Ep.ge ...
               blank ...
               blank];
else
    basicout = nan(In.ngroup, 10);
end

isfilled = isnan(basicin) & ~isnan(basicout);
isoutofrange = false(size(isfilled));
isoutofrange(:,7) = basicout(:,7) < 0 | basicout(:,7) > 1;
       
isinonly = isnan(basicout) & ~isnan(basicin);

dispval = basicin;
dispval(isfilled) = basicout(isfilled);

[nrow, ncol] = size(dispval);
    
gwidth = max(cellfun(@length, In.name));

dispstr = arrayfun(@(x) num2str(x,'%g'), dispval, 'uni', 0);
for ic = 1:ncol
    wdth(ic) = max(cellfun(@length, dispstr(:,ic)));
end
[dispstr{strcmp(dispstr, 'NaN')}] = deal(' ');

if cflag
    fmt = sprintf('%%-%ds ', 3, gwidth, wdth);
    cprintf('-text', fmt, '   ', 'Name', 'TL', 'HA', 'BH', 'B', 'PB', 'QB', 'EE', 'GE', 'GS', 'DI');
    fprintf('\n');

    for ir = 1:nrow
        fmt = sprintf('%%2d: %%-%ds', gwidth+1);
        fprintf(fmt, ir, In.name{ir});
        for ic = 1:ncol
            fmt = sprintf('%%%ds', wdth(ic)+1);
            if isoutofrange(ir,ic)
                cprintf('red', fmt, dispstr{ir,ic});
            elseif isfilled(ir,ic)
                cprintf('blue', fmt, dispstr{ir,ic});
            else
                cprintf('text',fmt, dispstr{ir,ic});
            end
        end
        fprintf('\n');
    end
else
    fmt = sprintf('%%-%ds ', 3, gwidth, wdth);
    fprintf(fmt, '   ', 'Name', 'TL', 'HA', 'BH', 'B', 'PB', 'QB', 'EE', 'GE', 'GS', 'DI');
    fprintf('\n');

    for ir = 1:nrow
        fmt = sprintf('%%2d: %%-%ds', gwidth+1);
        fprintf(fmt, ir, In.name{ir});
        for ic = 1:ncol
            fmt = sprintf('%%%ds', wdth(ic)+1);
            if isoutofrange(ir,ic)
                fprintf(fmt, dispstr{ir,ic});
            elseif isfilled(ir,ic)
                fprintf(fmt, dispstr{ir,ic});
            else
                fprintf(fmt, dispstr{ir,ic});
            end
        end
        fprintf('\n');
    end
end

           
           


