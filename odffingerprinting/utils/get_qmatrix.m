%% Read a q-matrix and calculate the RDSI reconstruction matrix

% Steven Baete
% NYU SOM CBI
% November 2016

%   Input:
%       basename: file with the diffusion directions for each shell
%       bvals: b-values of the shells
%       Delta: diffusion time, can be left empty
%       delta: gradient duration, can be left empty
%   Options (in 'keyword','value' -format)
%       ODFFile: tesselation file to use (default odf8.mat)
%       od_faces: faces of tesselation, to speed up execution
%       odf_vertices: vertices of tesselation, to speed up execution
%       Edge: edge factor for the RDSI reconstruction (default 1.25)
%   Output:
%       q: q-matrix
%       F: RDSI reconstruction matrix
%       odf_faces: faces of tesselation
%       odf_vertices: vertices of tesselation
%       qlist: q-values of the different shells
%       b_matrix: b-matrix

function [q,F,odf_faces,odf_vertices,qlist,b_matrix] = get_qmatrix(basename,bvals,Delta,delta,varargin)

if ((nargin < 3) | isempty(Delta)), Delta = 1; end;
if ((nargin < 4) | isempty(delta)), delta = 0; end;

odf_faces = [];
odf_vertices = [];
odffile = 'odf8.mat';
edge = 1.25;
if (nargin > 4)
    for i=1:2:size(varargin,2)
        if (~ischar(varargin{i}))
            error('get_qmatrix: Please use string-value pairs for input');
        end;
        switch varargin{i}
            case 'odf_faces'
                odf_faces = varargin{i+1};
            case 'odf_vertices'
                odf_vertices = varargin{i+1};
            case 'ODFfile'
                odffile = varargin{i+1};
            case 'Edge'
                edge = varargin{i+1};
            otherwise
                error(['get_qmatrix: I did not recognize the input-string ' varargin{i}]);
        end;
    end;
end;

% calculate q
qi = dlmread(basename,'\t');
if (size(qi,2) == 1)
    qi = dlmread(basename);
end;
qi = qi /max(qi(:));

nb = length(bvals);

if (length(Delta)==1), Delta = Delta*ones(nb,1);  end;
if (length(delta)==1), delta = delta*ones(nb,1);  end;

q = [];qDelta = [];qdelta = [];
for i=1:nb
    qlist(i) = sqrt(bvals(i)/(Delta(i)-delta(i)/3));
    q=vertcat(q,qi*qlist(i));   
    qDelta = vertcat(qDelta,ones(size(qi,1),1)*Delta(i));
    qdelta = vertcat(qdelta,ones(size(qi,1),1)*delta(i));
end

q=vertcat([0 0 0],q);
qDelta = vertcat(qDelta(1),qDelta);
qdelta = vertcat(0,qdelta);

% figure;plot3(q(:,1),q(:,2),q(:,3),'o');

% calculate F
if (isempty(odf_faces) | isempty(odf_vertices))
    load(odffile);
end;
b_matrix = cat(2,(sum(q.^2,2).*(qDelta-qdelta/3)),q./repmat(sqrt(sum(q.^2,2)), [1 3]))';
b_matrix(isnan(b_matrix)) = 0.0;

[F] = radial_recon_matrix(odf_vertices, b_matrix, ...
        edge,'Delta',qDelta,'delta',qdelta);
