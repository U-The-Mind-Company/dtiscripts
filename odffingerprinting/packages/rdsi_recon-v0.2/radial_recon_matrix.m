function [F E] = radial_recon_matrix(odfvertices, btable, edgefactor, varargin)
%RADIAL_RECON_MATRIX Create an encoding matrix to transform q-space image
%  data to an ODF.
%  Usage:
%    data2odf = radial_recon_matrix(...);
%    odf = data2odf * data;
%INPUTS
%  odfvertices - 3xNO sampling points for evaluating the ODF.  The 3 rows
%    of this matrix correspond to x, y, and z coordinates and should have a
%    max magnitude of 1.0.
%  btable - 4xNB matrix describing the q-space sampling scheme.  The 1st
%    row is b-value (s/mm^2) and rows 2-4 are the x, y, and z coordinates
%    of unit vectors giving the sampling orientation.
%
%OUTPUTS
%  F - NOxNB matrix that transforms q-space image intensities to ODF
%    samples.
%
%F Boada, S Yutzy
%University of Pittsburgh
%2/2/2012
% 
% Adapted by Steven Baete, NYU SOM CBI, 2013-2017

%% INPUT VALIDATION

if (nargin < 1 || isempty(odfvertices))
   error('odfvertices is required');
end;
if (ndims(odfvertices) ~= 2 || (size(odfvertices,1) ~= 3 && size(odfvertices,1) ~= 2))
   error('odfvertices must be a 3xNO matrix');
end;
if (nargin < 2 || isempty(btable))
   error('btable is required');
end;
if (ndims(btable) ~= 2 || (size(btable,1) ~= 4 && size(btable,1) ~= 3))
   error('btable must be a 4xNB matrix');
end;
if (nargin < 4 || isempty(edgefactor))
   edgefactor = 1.0;
end;
if (~isscalar(edgefactor))
   error('edgefactor must be a scalar');
end;

Delta = 1;
delta = 0;
bcomplete = [];
if (nargin > 4)
    for i=1:2:size(varargin,2)
        if (~ischar(varargin{i}))
            error('radial_recon_matrix: Please use string-value pairs for input');
        end;
        switch varargin{i}
            case 'bcomplete'
                bcomplete = varargin{i+1};        
            case 'Delta'
                Delta = varargin{i+1};
            case 'delta'
                delta = varargin{i+1};
            otherwise
                error(['radial_recon_matrix: I did not recognize the input-string ' varargin{i}]);
        end;
    end;
end;

dimi = size(odfvertices,1);

%processing is significantly faster if we transpose the matrices
%the only reason the matrices have the dimensions they do is because that
%is the size that DSI-Studio expects
odfvertices_T = odfvertices';
btable_T = btable';
if ~isempty(bcomplete)
    btable_TC = bcomplete';
    Q_TABLE_C = btable_TC(:,1+(1:dimi)) .* repmat(sqrt(btable_TC(:,1)/(Delta-delta/3)), [1 dimi]);
    dcfcomp = dcfcalc(Q_TABLE_C);
    diffs = diff(sqrt(sort((btable_TC(:,1))/(Delta-delta/3))));
    stepsC = floor(diffs/(max(diffs/20)))*(max(diffs/20));
    stepsC = sort(stepsC(stepsC > 0));
    qvalcompl = sqrt(sum(Q_TABLE_C.^2,2));
else
    dcfcomp = [];
end;

%% PROCESSING
%Convert to q-space representation
Q_TABLE = btable_T(:,1+(1:dimi)) .* repmat(sqrt(btable_T(:,1)/(Delta-delta/3)), [1 dimi]);

diffs = diff(sqrt(sort((btable_T(:,1))/(Delta-delta/3))));
steps = floor(diffs/(max(diffs/20)))*(max(diffs/20));
steps = sort(steps(steps > 0));

qval = sqrt(sum(Q_TABLE.^2,2));

% normalize the dcf
if (~isempty(dcfcomp))
    dcf = dcfcomp(1:size(Q_TABLE,1))/sum(dcfcomp)*sqrt((qvalcompl(end)/qval(end)));
    Q_STEP = stepsC(end);
else
    dcf = dcfcalc(Q_TABLE);	
    dcf = dcf/sum(dcf);
    Q_STEP = steps(end);
end;
h = dcf;

Rmax = edgefactor * 1/Q_STEP/2;

display([' radial_recon_matrix: nshell ' num2str(length(steps)) ' Rmax Q_step ' num2str(edgefactor * 1/Q_STEP/2)]);

E = odfvertices_T * Q_TABLE';

F = -sincpp(2*pi.*E * Rmax).*repmat((h)',[size(E,1),1]);

