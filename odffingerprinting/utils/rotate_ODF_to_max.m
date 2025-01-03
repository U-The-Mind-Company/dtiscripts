%% rotate the maximum point on the ODF to the Z-axis
%   invert the rotation if maxind and R supplied

% Steven Baete
% NYU SOM CBI
% October 2016

%   Input:
%       odf: odf-values
%       odf_faces, odf_vertices: odf tesselation
%       maxind, R: supply to invert the rotation
%       precise: more precise rotation than based on the tesselation
%   Output:
%       odfrot: values of the odf rotations
%       maxind: index of rotation matrix to use for each odf
%       R: rotation matrices used

function [odfrot,maxind,R] = rotate_ODF_to_max(odf,odf_vertices,maxind,R,precise)

if (nargin == 4) invert = true;
else invert = false; end;
if (nargin < 5 | isempty(precise))
    precise = true;
end;

% maximum values
if (~invert)
    [~,maxind] = max(odf,[],2);
end;

sodf = size(odf);

odf = cat(2,odf,odf);

%% filenames
randnum = num2str(round(rand(1)*1000));
shname = ['tmp' randnum '_sh.nii'];
peakname = ['tmp' randnum '_peak.nii'];

%% calculate the rotation matrices
if ~invert
    for i = 1:(size(odf_vertices,2)/2)
        R(i).R = rotation_matrix_twovectors(odf_vertices(:,i),[0 0 1]');
    end;
end;

%% convert odf to sh

shtrans = shtransformation(odf_vertices',14);
shorig = amp2sh(odf,shtrans);
clear odf;

if (precise)
    reshodf = save_to_mrtrix_nii(shorig,shname);
    system(['sh2peaks -quiet -force -num 1 ' ...
        shname ' ' peakname]);
    peaks = load_from_mrtrix_nii(peakname,reshodf);
    system(['rm ' shname ' ' peakname]);
    negz = (peaks(:,3)<0);
    peaks(negz,:) = -peaks(negz,:);
    peaks = normalizevector(peaks);
    maxind(:) = 0;
    l = 0;
    for i = 1:length(maxind)
        if (maxind(i) == 0)
            l = l+1;
            sel = (abs(angle_twovectors(peaks(i,:),peaks)) < 0.5) & maxind == 0;
            maxind(sel) = l;
            R(l).R = rotation_matrix_twovectors(peaks(i,:),[0 0 1]');
        end;
    end;
    display(['  rotate_ODF_to_max needs ' num2str(l) ' rotations']);
end;

%% rotation
odfrot = zeros(sodf);

for i = 1:length(R)
    % which odfs are we talking about
    ind = (maxind==i);
    
    if (sum(ind) > 0)        
        % rotation matrix
        if invert
            anglesrot = R(i).R*odf_vertices;
        else
            anglesrot = R(i).R'*odf_vertices;
        end;
        
        shtrans = shtransformation(anglesrot',14);
        odft = sh2amp(shorig(ind,:),shtrans);
        odfrot(ind,:) = odft(:,1:sodf(2));
    end;
end;

% remove negative values
odfrot(odfrot(:) < 0) = 0;
