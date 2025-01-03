%% Find the peaks of ODFs using local-maximum search
%  based on the find_ODF_peak-script by Frank Yeh, DSIStudio

% Steven Baete
% NYU SOM CBI
% November 2016

%   Input:
%       odf: odf-values
%       odf_faces, odf_vertices: odf tesselation
%       MAX_FIBERS: maximum number of fibers to search for (default 5)
%   Output:
%       dirs: identified fiber directions
%       qa: identified assosciated qa-values

function [dirs,qa] = find_ODF_peak(odf,odf_faces,odf_vertices,MAX_FIBERS)

if(nargin < 4 | isempty(MAX_FIBERS))
    MAX_FIBERS = 5;
end;

[nodf,nvert] = size(odf);
dirs = NaN*ones(nodf,MAX_FIBERS,2);
qa = zeros(nodf,MAX_FIBERS);

odf_faces = odf_faces + 1;
odf_faces = odf_faces - (odf_faces > (2*nvert))*(2*nvert);

angles = odf_vertices(:,1:nvert);
[azimuth,elevation,~] = cart2sph(angles(1,:),angles(2,:),angles(3,:));

for i = 1:nodf
    % the odf
    odft = [odf(i,:) odf(i,:)];
    is_peak = odft;
    temp1 = odft(odf_faces(1,:));
    temp2 = odft(odf_faces(2,:));
    temp3 = odft(odf_faces(3,:));
    temp = temp2 > temp1 | temp3 > temp1 | (temp2 == temp1 & temp3 == temp1);
    is_peak(odf_faces(1,temp)) = 0;
    temp = temp1 > temp2 | temp3 > temp2 | (temp1 == temp2 & temp3 == temp2);
    is_peak(odf_faces(2,temp)) = 0;
    temp = temp2 >= temp3 | temp1 >= temp3;
    is_peak(odf_faces(3,temp)) = 0;
    [values,ordering] = sort(-is_peak);
    p = ordering(values < 0);
    
    % clean and sort the peaks
    pks = p(1:2:end);
    np = length(pks);
    if (np > MAX_FIBERS)
        np = MAX_FIBERS;
        pks = pks(1:np);
    end;
    pks(pks > nvert) = pks(pks > nvert) - nvert;
    
    %calculate qa
    odfVals = odft(pks);
    %sort qa by odf value to arrange into fa0, fa1, etc.
    [odfVals, idx] = sort(odfVals,'descend');
    pks = pks(idx);
    qat = odfVals - min(odft);
    
    % save results
    qa(i,1:np) = qat;
    dirs(i,1:np,1) = azimuth(pks);    
    dirs(i,1:np,2) = elevation(pks);
end;