%% save dat to filename for processing with MRtrix3

% Steven Baete
% NYU SOM CBI
% October 2016

function resh = save_to_mrtrix_nii(dat,filename,header,datatype)

if (nargin < 3),     header = []; end;
if (nargin < 4),     datatype = 64; end;

maxsize = 30000;
so = size(dat);
added = 0;
if (length(so) == 2)
  mf = factor(so(1));
  while(max(mf) > maxsize) % avoid any dimension larger than 2^15
      dat = cat(1,dat,ones(1,so(2)));
      mf = factor(size(dat,1));
      added = added + 1;
  end
  m = zeros(size(dat,1),1,1,size(dat,2));
  m(:,1,1,:) = dat;
  ms = [1,1,1];
  mo = size(m);ms(1:length(mo)) = mo;
  mf = factor(ms(1));
  ind = find(cumprod(mf)< maxsize,1,'last');
  if (isempty(ind)) ind = length(mf); end;
  ms(1) = prod(mf(1:ind));
  if (ind < length(mf))
      ind2 = find(cumprod(mf((ind+1):end))< maxsize,1,'last');
      if (isempty(ind2)) ind2 = length(mf)-ind; end;
      ms(2) = prod(mf(ind+(1:ind2)));
      if (ind2+ind < length(mf))
          ind3 = find(cumprod(mf((ind+ind2+1):end))< maxsize,1,'last');
          if (isempty(ind3)) ind3 = length(mf)-ind2-ind; end;
          ms(3) = prod(mf(ind+ind2+(1:ind3)));
          if (ind3+ind2+ind < length(mf))
              error;
          end;
      end;
  end;
  if (ms(2)==1 & length(factor(ms(1))>1))
      mf = factor(ms(1));
      [ms(1),ind] = max(mf);
      ms(2) = prod(mf([1:(ind-1),(ind+1):end]));
  end;
  if (ms(3)==1 & length(factor(ms(2))>1))
      mf = factor(ms(2));
      [ms(2),ind] = max(mf);
      ms(3) = prod(mf([1:(ind-1),(ind+1):end]));
  end;
  dat = reshape(m,ms);
  resh = so;
  sq = true;
else
  resh = [];
  sq = false;
end;

if (~isempty(header))
    header.img = dat;
    header.hdr.dime.dim(1) = 4;
    header.hdr.dime.dim(5) = so(4);
    save_untouch_nii(header,filename);
else
    nii = make_nii(dat,[],[],datatype);
    if (length(size(dat)) == 2)
        nii.hdr.dime.dim(1) = 3;
    end;
    save_nii(nii,filename);
end;


