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
odfname = ['tmp' randnum '_odf.nii'];
shname = ['tmp' randnum '_sh.nii'];
shpname = ['tmp' randnum '_shp_XX.nii'];
odfrotname = ['tmp' randnum '_odfrot_XX.nii'];
dirname = ['tmp' randnum '_XX.dirs'];
peakname = ['tmp' randnum '_peak.nii'];

%% write directions to file
% so that it agrees with mrtrix
angles = odf_vertices;
[azimuth,elevation,~] = cart2sph(angles(1,:),angles(2,:),angles(3,:));
fid = fopen(dirname,'w');
fprintf(fid,'%f %f\n',[azimuth;pi/2-elevation]);
fclose(fid);

%% calculate the rotation matrices
if ~invert
    for i = 1:(size(odf_vertices,2)/2)
        R(i).R = rotation_matrix_twovectors(odf_vertices(:,i),[0 0 1]');
        R(i).fname = strrep(dirname,'XX',[num2str(i,'%03.0f')]);
    end;
end;

%% convert odf to sh
% mrtrix amp2sh
reshodf = save_to_mrtrix_nii(odf,odfname);
% clear odf;

system(['amp2sh -quiet -force -lmax 14 ' ...
    '-directions ' dirname ...
    ' ' odfname ' ' shname]);

if (precise)
    system(['sh2peaks -quiet -force -num 1 ' ...
        shname ' ' peakname]);
    peaks = load_from_mrtrix_nii(peakname,reshodf);
    system(['rm ' peakname]);
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
            R(l).fname = strrep(dirname,'XX',[num2str(i,'%03.0f')]);
        end;
    end;
    display(['  rotate_ODF_to_max needs ' num2str(l) ' rotations']);
end;

%% load the sh-coefficients
shorig = load_from_mrtrix_nii(shname,reshodf);

%% rotation
odfrot = zeros(sodf);

comstr = [];
l = 1;
for i = 1:length(R)
    % which odfs are we talking about
    ind = (maxind==i);
    
    if (sum(ind) > 0)        
        % write sh to file
        reshsh(i).re = save_to_mrtrix_nii(shorig(ind,:),strrep(shpname,'XX',num2str(i)));
        
        % rotation matrix
        if invert
            anglesrot = R(i).R*odf_vertices;
        else
            anglesrot = R(i).R'*odf_vertices;
        end;
        [azimuth,elevation,~] = cart2sph(anglesrot(1,:),anglesrot(2,:),anglesrot(3,:));
        fid = fopen(R(i).fname,'w');
        fprintf(fid,'%f %f\n',[azimuth;pi/2-elevation]);
        fclose(fid);
                
        % matrix sh2amp
        comstr = [comstr 'sh2amp -quiet -force ' ...
            ' ' strrep(shpname,'XX',num2str(i)) ...
            ' ' R(i).fname ...
            ' ' strrep(odfrotname,'XX',num2str(i)) ' \n '];%
    end;
    % --------------------------------
    if (length(comstr) > 100000*100 | i == length(R))
        fid = fopen('command.txt','wt');
        fprintf(fid,comstr);
        fclose(fid);
        system(['./utils/par_exec.sh command.txt shamp']);
        comstr = [];
    end;
    % --------------------------------
end;

% read tmp-files first
for i = 1:length(R)
    % which odfs are we talking about
    ind = (maxind==i);
    
    if (sum(ind) > 0)   
        % read the rotated odfs
        try
            tmp = load_from_mrtrix_nii(strrep(odfrotname,'XX',num2str(i)),reshsh(i).re);
            if (sum(tmp(:))==0)
                error
            end;
        catch
            error
        end;
        odfrot(ind,:) = tmp(:,1:sodf(2));
    end;
end;

% remove negative values
odfrot(odfrot(:) < 0) = 0;

%% clear files
system(['rm ' strrep(odfrotname,'XX','*') ' ' ...
            strrep(shpname,'XX','*') ' ' shname ...
    ' ' odfname ' ' strrep(dirname,'XX','*') ' command.txt']);
