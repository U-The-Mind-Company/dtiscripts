%% rotate directions according to R(maxind)
%   invert the rotation if invert is true

% Steven Baete
% NYU SOM CBI
% November 2016

function [dirrot] = rotate_dir(dir,maxind,R,invert)

if (nargin < 4)
    invert = false;
end;

dirrot = zeros(size(dir));
for i = 1:length(R)
    if (mod(i,10000) == 0)
        display(['    rotating dir : [' num2str(i) '/' num2str(length(R)) ']']);
    end;
    % which odfs are we talking about
    ind = (maxind==i);
    
%     if (i == 1) %identity matrix
%         dirrot(ind,:,:) = dir(ind,:,:);
%     else
        if (sum(ind) > 0)    
            % the directions to rotate
            dirs = dir(ind,:,:);
            s = size(dirs);
            dirslist = reshape(dirs,[s(1)*s(2),s(3)]);
            nzeroind = sum(isnan(dirslist),2)==0;
            dirslistshort = dirslist(nzeroind,:);

            xyz = zeros(size(dirslistshort,1),3);

            [xyz(:,1),xyz(:,2),xyz(:,3)]=sph2cart(dirslistshort(:,1),...
                dirslistshort(:,2),ones(size(dirslistshort,1),1));

            % rotation
            if invert
                xyzrot = R(i).R'*xyz';
            else
                xyzrot = R(i).R*xyz';
            end;

            dirslistrotshort = zeros(size(dirslistshort));
            [dirslistrotshort(:,1),dirslistrotshort(:,2),~] = ...
                cart2sph(xyzrot(1,:)',xyzrot(2,:)',xyzrot(3,:)');

            dirslistrot = NaN*zeros(size(dirslist));
            dirslistrot(nzeroind,:) = dirslistrotshort;
            dirsrot = reshape(dirslistrot,[s(1),s(2),s(3)]);

            dirrot(ind,:,:) = dirsrot;
        end;
%     end;
end;