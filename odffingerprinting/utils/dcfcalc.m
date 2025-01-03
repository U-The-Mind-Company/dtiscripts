% Steven Baete
% NYU SOM CBI
% September 2017
% density compensation function, this extends h calculated above
%    the dcf is a multiplication of 3 dimensions
%       two dimensions (orthogonal to the radial lines) are each similar to the
%       2D-radial dcf, hence, ~ (qval)^2
%       the third dimension, along the radial lines, depends on the spacing
%       of the samples.

% INPUTS
%   qtable - Nx3 matrix with the q-space sampling scheme, each row contains
%            one sampling orientation, scaled by q-value 
%   estmethod - 'theor' (default),'geom','pipe2','pipe','voronoi','jackson': method
%            to estimate the dcf. All but 'theor' and 'geom' need the
%            nufft-package (http://web.eecs.umich.edu/~fessler/irt/irt/nufft) in the path.
%         'theor' and 'geom' are variations of theoretical volumetric
%            calculations of the dcf, based on an RDSI q-space sampling
%            (i.e. multiple samples on each of a number of radially outward
%            lines in q-space)
%            i.e. The dcf is a multiplication of 3 dimensions; two dimensions 
%               (orthogonal to the radial lines) are each similar to the  2D-radial dcf, 
%               hence, ~ (qval)^2; the third dimension, along the radial lines, 
%               depends on the spacing of the samples.
%         'pipe2','pipe','voronoi','jackson' are variations of methods used
%            for non-uniform fft dcf calculation. These are more general
%            regarding the assumptions of the q-space sampling. See
%            documentation in the nufft-package
% OUTPUTS
%   dcf -  Nx1 matrix containing the density compensation function

function dcf = dcfcalc(qtable,estmethod)

if (nargin < 2 | isempty(estmethod))
    estmethod = 'theor';
end;
normalize = true;

dimi = size(qtable,2);

qval = sqrt(sum(qtable.^2,2));
[qshells,nmeasshell,shell] = bval2shells(qval');
steps = qshells;

% k-space traj and FoV (for normalization)
nz = 2;
dq = mean(diff(qshells*1e3)); % m-1
Rmax = 1/dq; % m
FoV = 2*Rmax*[1,1,1];
res = (2*length(qshells)+1)*[1,1,1];
kspace = qtable*FoV(1)*1/res(1)*(pi/2);
kspace = [kspace;-kspace];
kspace = kspace/max(kspace(:))*1/2*2*pi;

% Nufft
if (~exist('nufft_init','file'))
    if (~strcmp(estmethod,'theor') & ~strcmp(estmethod,'geom'))
        warning(['The IRT package seems not installed. You will find it at:' ...
            ' http://web.eecs.umich.edu/~fessler/irt/irt/nufft .' ...
            ' dcfcalc.m will use a simple calculation of the DCF.']);
        estmethod = 'theor';
    end;
    normalize = false;
    warning(['The IRT package seems not installed. dcfcalc.m will not normalize the DCF.']);
else
    st = nufft_init(kspace,...
        res, [5,5,5], 2*res,res, 'minmax:kb');
end;

% calculate dcf
switch estmethod
    case 'theor'
        % the dcf is a multiplication of 3 dimensions
        % two dimensions (orthogonal to the radial lines) are each similar to the
        % 2D-radial dcf, hence, ~ (qval)^2
        % the third dimension, along the radial lines, depends on the spacing
        % of the samples.
        dcf = qval/max(qval);
        if (dimi == 3)
            dcf = dcf.*dcf;
        end;
        dd = diff(steps)';
        dd=dd/sum([dd;dd(end)/2])*(length(dd)+1);
        dcfw = 1/2*[0;dd]+1/2*[dd;dd(end)];
        for i = 1:length(steps)
            sel = (shell == i);
            dcf(sel) = dcf(sel)*dcfw(i);
        end;
    case 'geom'
        % volume in q-space, or, (4/3 pi r2^3 - 4/3 pi r1^3) / nmeas with r1
        % and r2 dividing the space between the neighbouring shells
        diffsteps = diff(steps)';
        dsplit = [0,(qshells(1:(end-1))+qshells(2:end))/2,qshells(end)+diffsteps(end)/2];
        dsplit(2) = dsplit(2)/10;
        dcfshells = 4/3*pi*(dsplit(2:end).^3 - dsplit(1:(end-1)).^3)./nmeasshell;
        dcfshells = dcfshells/dcfshells(end);
        dcf = dcfshells(shell)';
        
    % below are for non-uniform direction distribution (general case)
    case 'pipe2'
        % add extra outside shell
        qvalextra = (length(qshells)*qshells(2:end)/(1:(length(qshells)-1)));
        dirextra = dlmread('dir1000uniform.txt');
        dirextra = qtable((shell == length(qshells)),:);
        kspacet = [kspace;dirextra*1/2*2*pi*qvalextra/qshells(end)]*qshells(end)/qvalextra;
        st3 = nufft_init(kspacet/nz,...
            nz*res, [5,5,5], 2*nz*res,res, 'minmax:kb');
        H.arg.st = st3;
        Dest = ir_mri_density_comp_v2(kspacet, estmethod,...
              'G',H,'arg_pipe',{'Nsize',nz*res(1),'niter',50});
        Dest = Dest(1:(end-size(dirextra,1)));
        dcf = Dest(1:(end/2));
%     case 'pipe2'
%         st2 = nufft_init(kspace/nz,...
%             nz*res, [5,5,5], 2*nz*res,res, 'minmax:kb');
%         G.arg.st = st2;
%         Dest = ir_mri_density_comp_v2(kspace, estmethod,...
%               'G',G,'arg_pipe',{'Nsize',nz*res(1),'niter',50});
%         dcf = Dest(1:(end/2));
    case 'voronoi'
        % add extra outside shell
        qvalextra = (length(qshells)*qshells(2:end)/(1:(length(qshells)-1)));
        dirextra = dlmread('dir1000uniform.txt');
        kspacet = [kspace;dirextra*1/2*2*pi*qvalextra/qshells(end)]*qshells(end)/qvalextra;
        st3 = nufft_init(kspacet/nz,...
            nz*res, [5,5,5], 2*nz*res,res, 'minmax:kb');
        H.arg.st = st3;
        Dest = ir_mri_density_comp_v2(kspacet, estmethod,...
          'G',H,'fix_edge',0);
        Dest = Dest(1:(end-size(dirextra,1)));
        dcf = Dest(1:(end/2));
    case {'jackson';'pipe'}
        st2 = nufft_init(kspace/nz,...
            nz*res, [5,5,5], 2*nz*res,res, 'minmax:kb');
        G.arg.st = st2;
        Dest = ir_mri_density_comp_v2(kspace, estmethod,'G',G);
        dcf = Dest(1:(end/2));
end;

% normalize 
if normalize
    out = (st.p'*(repmat(dcf,[2,1]).*((st.p)*ones(prod(2*res),1))));
    dcf = dcf/mean(abs(out(:)));
end;
