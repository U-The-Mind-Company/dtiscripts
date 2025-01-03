%% load dat from filename for processing with MRtrix3

% Steven Baete
% NYU SOM CBI
% October 2016

function dat = load_from_mrtrix_nii(filename,resh)

if (nargin < 2), resh = []; end;

nii = load_nii(filename);
dat = double(nii.img);

if (~isempty(resh))
    s = size(dat);
    resh(end) = s(end);
    if(prod(s) == prod(resh(1:(end-1))))
        resh(end) = 1;
    end;
    try
        dat = reshape(dat,resh);
    catch
        dat = reshape(dat,[prod(s(1:3)),resh(2)]);
        dat = dat(1:resh(1),:);
    end;
end;