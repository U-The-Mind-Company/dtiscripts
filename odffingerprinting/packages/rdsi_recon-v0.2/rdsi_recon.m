function [fibfileout] = rdsi_recon(srcfile,varargin)
%rdsi_recon Perform radial q-space sampling reconstruction on a
%  DSI-Studio format image .src file
% INPUTS
%  srcfile - Path to the .src or .src.gz file created by DSI-Studio. If the
%    file specified is in gzipped format, it will be decompressed to a temporary
%    .src file, which will be removed at the end of the function.
% OPTIONAL INPUTS (in string-value pairs, e.g. " 'maxfibers',5, " )
%  fibfile - Path in which to save the .fib.gz output file created by this
%    function.  If this is not specified (or empty), defaults to the filename of
%    srcfile with '.fib.gz' as extension.
%  maskfile - Path to a DSI-Studio format mask file.  If this is not specified
%    (or empty), the script will generate a mask and save it as 
%       srcfile.mask.txt .
%  odftess - Name of the ODF sampling file to use.  Valid values are:
%    'odf4' - 162 vertices
%    'odf5' - 252 vertices
%    'odf6' - 362 vertices
%    'odf8' - 642 vertices (default)
%  maxfibers - Maximum number of fibers to record in the output file.
%    maxfibers must be >= 1, and defaults to 5;
%  output - If true, progress messages are output to the screen throughout the
%    processing.  Defaults to true.
%  edgefactor - Factor used to adjust the maximum displacement distance Rm
%    to allow for ODF sharpening. A value higher than 1 will sharpen the ODF,
%    lower than 1 will blunt the ODF, default 0.8.
%  phantommask - to use simple thresholding instead of brain segmenation
%    for masking
%  flipBaxis - flips the listed axis of the b_table
%OUTPUTS
%  fibfile - Name of the .fib.gz file written
%
% -----------------------------------------------------------------------
% This method has been published in:
%  Baete, Yutzy and Boada. Radial q-space sampling for DSI.
%       Magn Reson Med, 76(3):769-80, 2016
%  http://onlinelibrary.wiley.com/doi/10.1002/mrm.25917/abstract
%
%  and expanded on in:
%  Baete and Boada. Accelerated Radial Diffusion Spectrum Imaging using 
%       a Multi-Echo Stimulated Echo diffusion sequence.
%	Magn. Reson. Med., 79(1):306-316, 2018.
%  http://onlinelibrary.wiley.com/doi/10.1002/mrm.26682/abstract
% -----------------------------------------------------------------------
%
% Steven Baete
% NYU SOM CBI
% 02/23/2015, adapted September 2017
%
% original script from: Stephen Yutzy, University of Pittsburgh
% 02/15/2012

%% PARAMETER CHECKING
if (~exist(srcfile, 'file'))
    error('rdsi_recon: Cannot find srcfile: %s', srcfile);
end

fibfile = [srcfile '.fib'];
maskfile = [srcfile '.data.txt'];
odftess = 'odf8';
maxfibers = 5;
edgefactor = 0.8;
output = true;
phantommask = false;
flipbaxis = [];
if (nargin > 1)
    for i=1:2:size(varargin,2)
        if (~ischar(varargin{i}))
            error('rdsi_recon: Please use string-value pairs for input');
        end
        switch varargin{i}
            case 'fibfile'
                fibfile = varargin{i+1};
            case 'odftess'
                odftess = varargin{i+1};
            case 'maskfile'
                maskfile = varargin{i+1};
            case 'maxfibers'
                maxfibers = varargin{i+1};
            case 'output'
                output = varargin{i+1};
            case 'edgefactor'
                edgefactor = varargin{i+1};
            case 'phantommask'
                phantommask = varargin{i+1};
            case 'flipBaxis'
                flipbaxis = varargin{i+1};
            otherwise
                error(['rdsi_recon: I did not recognize the input-string ' varargin{i}]);
        end
    end
end

%% SET UP VARIABLES
[~,src_filename,src_fileext] = fileparts(srcfile);
if (strcmpi(src_fileext, '.gz'))
   if (output);fprintf(1, 'Decompressing src file %s%s...', src_filename, src_fileext);end
   srcfileorig = srcfile;
   if (isunix())
       unix(['gunzip -c ' srcfileorig ' > ' srcfileorig(1:(end-3))]);
   else
       gunzip(srcfileorig);
   end
   srcfile = srcfileorig(1:(end-3));
   if (output);fprintf(1, 'done\n');end
end

if (output);fprintf(1, 'Setting up...');end
%load dimension, voxel_size, b_table.  Images will be loaded one at a time
%to save memory. Testing shows that this increases the runtime of the
%script by about 10%, which is an acceptable tradeoff for cutting the
%memory footprint by 2x.
load(srcfile, '-mat', 'b_table','dimension', 'voxel_size');

% different convention
if (~isempty(flipbaxis))
    for i = 1:length(flipbaxis)
        b_table(flipbaxis(i)+1,:) = -b_table(flipbaxis(i)+1,:);
    end;
end;

odfVarSize = 20000; %maximum size of an ODF variable in the output file

%load odf_vertices and odf_faces for our chosen odf sampling resolution
load(odftess);

odfHalfLength = size(odf_vertices,2)/2;
% calculate this here, so that it doesn't have to be calculated time and
% time again in find_peak_sb
odf_faces_fp = odf_faces + 1;
odf_faces_fp = odf_faces_fp - (odf_faces_fp > 2*odfHalfLength)*(2*odfHalfLength);

if (output);fprintf(1, 'done\n');end

%% PROCESS ODFS
%load and normalize images
if (output);fprintf(1, 'Loading and normalizing images from %s...', src_filename);end
%process echo normalization
idxb0 = find(b_table(1,:) < 125); 
%check for b=0 images outside of the beginning of the b-table
if (length(idxb0) > 1)
   warning('found multiple b=0 images'); 
end

%load images
img = zeros([size(b_table,2) dimension]);
for qq = 1:size(b_table,2)
   %determine the name of the imageN variable
   imgName = sprintf('image%d', qq-1);
   load(srcfile, '-mat', imgName);
   eval(sprintf('img(qq,:,:,:) = reshape(%s, dimension);', imgName));
   clear(imgName);
end

%% calculate mask and csfmask

img = permute(img,[2,3,4,1]);
img(isnan(img(:))) = 0;

[mask,~,masks] = segmentbrain(img,maskfile,0,[],...
            floor(size(img,3)*2/5):floor(3/5*size(img,3)));
if (phantommask)
    maskc = img(:,:,:,1) > 200;
    [xt,yt,zt] = meshgrid(1:dimension(2),1:dimension(1),1:dimension(3));
    maskt = reshape(maskc,[1,prod(dimension)]);
    mask = [yt(maskt == 1);xt(maskt == 1);zt(maskt == 1)]';
    masks{1} = maskc;masks{2}=maskc;masks{3}=maskc;
end;
img = permute(img,[4,1,2,3]);
csfmask = double(masks{1});
% which points of the mask are csf?
csfinmask = csfmask(sub2ind(dimension,mask(:,1),mask(:,2),mask(:,3)));

% mask3d = zeros(dimension);
% mask3d(sub2ind(dimension,mask(:,1),mask(:,2),mask(:,3))) = 1;
% figure;imagesc(montageSB(mask3d))

%% normalize images
        
% save the b0 for later
img0 = permute(mean(img(idxb0,:,:,:),1),[2,3,4,1]);

% remove extra b0-images
if (length(idxb0) > 0)
    img(idxb0(1),:,:,:) = img0;
    img(idxb0(2:end),:,:,:) = [];
    b_table(:,idxb0(1)) = mean(b_table(:,idxb0),2);
    b_table(:,idxb0(2:end)) = [];
    idxb0 = idxb0(1);
end;

% divide out the b0
for qq = 1:size(b_table,2)
    img(qq,:,:,:) = img(qq,:,:,:) ./ permute(img0,[4,1,2,3]);
    b_table(:,qq) = subtractDiffDir(b_table(:,qq),b_table(:,idxb0));
end
img((img>1)) = 1;

% weigh with b0
img = img.*repmat(permute(img0,[4,1,2,3]),[size(img,1),1,1,1]);
b0 = reshape(img0, [1 prod(dimension)]);

%% find the shells
table = b_table(1,:);
[sbtable,~] = sort(table);
inds = find(diff(sbtable)./max(sbtable(1:(end-1)),50) > 0.1);
shellstep = (sbtable(inds)+sbtable(inds+1))/2;
bvalues = sbtable(inds+1);
nmeasshell = diff([inds,length(sbtable)]);
nshells = length(shellstep);
shell = zeros(size(sbtable));
indc = 1:size(b_table,2);
for i = 1:nshells
    shell(table > shellstep(i)) = i;
end
for i = 1:nshells
    shellind{i} = indc(shell <= i);
end
for i = 0:nshells
    pershellind{i+1} = indc(shell == i);
    pershellindn(i+1) = length(pershellind{i+1});
end

if (output);fprintf(1, 'done\n');end

%% reconstruct ODFs

Ipsi0 = 0;
%create data -> ODF encoding matrix
if (output);fprintf(1, 'Creating encoding matrix...\n');end
% normalize the b_table vectors input for the
% radial_recon_matrix-calculation
b_table2 = b_table;
b_table2(2:4,:) = b_table2(2:4,:)./ repmat(sqrt(sum(b_table2(2:4,:).^2)),[3,1]);
b_table2(isnan(b_table2)) = 0;
for i = 1:nshells
    data2odf{i} = radial_recon_matrix(odf_vertices, b_table2(:,shellind{i}), ...
        edgefactor, 'bcomplete',b_table2);
    data2odf{i}=data2odf{i}/norm(data2odf{i});
end
if (output);fprintf(1, 'done\n');end


%reconstruct ODFs
if (output);fprintf(1, 'Reconstructing ODFs...');end
odf = zeros(odfHalfLength, size(mask,1));
shia = zeros(1, size(mask,1));
for mm = 1:size(mask,1)
    mt = mask(mm,:);
    data = img(:, mt(1), mt(2), mt(3));

    % average value per shell
    shi = nshells;
    for i = 0:nshells
        shellval(i+1) = mean(data(pershellind{i+1}));
    end            
    % how much shells are useful data? -> only use the useful data
    % a.k.a. where is the noise floor for this voxel?
    dshell = -diff(shellval);
    shi = max(find(dshell(2:end) < 0.1*dshell(1:(end-1)),1,'first'),2);
    if (isempty(shi))
        shi = nshells;
    end
    odfTmp = data2odf{shi} * data(shellind{shi});

   odfTmp(isnan(odfTmp)) = 0;
   mmin = min(odfTmp);
   if (mmin < 0 ) %$% changed position
       odfTmp((odfTmp <= 0)) = 0;
       mmin = 0;
   end       
   if (sum(odfTmp) == 0) odfTmp(:) = 1e-8; end

   odf(:,mm) = odfTmp(1:odfHalfLength);
   shia(mm) = shi;
   % Ipsi0: see Frank Yeh, Generalized q-Space Sampling, IEEE TMI (29) 1626-35
   if (csfinmask(mm)) % this is csf
       minIpsi0(mm) = min(odfTmp);
       Ipsi0 = max(Ipsi0,min(odfTmp));
   end
end
if (output);fprintf(1, 'done\n');end

%find ODF peaks
if (output);fprintf(1, 'Finding ODF peaks...');end
peaks = cell(1, size(mask,1));
maxPeaks = 0;
for mm = 1:size(odf,2)
   peaksTmp = find_peak_sb(vertcat(odf(:,mm), odf(:,mm)), odf_faces_fp);
   if (length(peaksTmp)/2 > maxPeaks)
      maxPeaks = length(peaksTmp)/2;
   end
   peaks{mm} = peaksTmp;
end
if (output);fprintf(1, 'done\n');end

% plot diffusion directions
% figure;sel = (b_table(1,:)>3000);plot3(b_table(2,sel),b_table(3,sel),b_table(4,sel),'o');axis equal;xlabel('x');ylabel('y');zlabel('z');
% figure;sel = (b_table(1,:)>3000);hold on;plot(b_table(2,sel),ones(1,sum(sel)),'or');plot(b_table(3,sel),1.2*ones(1,sum(sel)),'og');plot(b_table(4,sel),1.4*ones(1,sum(sel)),'ob');ylim([0.5 1.5])

% figure;hold on;sel = peaksTmp;plot3(odf_vertices(1,:).*odf([1:end,1:end],mm)',odf_vertices(2,:).*odf([1:end,1:end],mm)',odf_vertices(3,:).*odf([1:end,1:end],mm)','ob');plot3(odf_vertices(1,sel),odf_vertices(2,sel),odf_vertices(3,sel),'or');axis equal;xlabel('x');ylabel('y');zlabel('z');

% scale z0, so that the ODF of free water (CSF here) is 1
% see Frank Yeh, Generalized q-Space Sampling, IEEE TMI (29) 1626-35
% z0 here is the inverse of Frank Yeh's z0
if (exist('minIpsi0'))
    Ipsi0 = mean(minIpsi0(csfinmask==1));
else
    Ipsi0 = 1;
end
vol = 1;
z0 = vol/Ipsi0;
if (z0 <= 0)
    display(['    rdsi_recon: Warning: Found z0 smaller or equal than 0 (' num2str(z0) '), put it to 1!']);
    z0 = 1;
end

%% WRITE RESULTS
if (output);fprintf(1, 'Saving fa and index output...');end
save(fibfile, 'dimension', 'odf_vertices', 'odf_faces', 'voxel_size','mask','csfinmask', '-mat', '-v4');

%create fa* and index* output variables based on maxfibers
for ff = 0:maxfibers-1
   fal{ff+1} = zeros(dimension);
   indexl{ff+1} = zeros(dimension);
end

%write fa and index to the appropriate variables
for mm = 1:size(mask,1)
    mt = mask(mm,:);
   %get the peaks and odf values
   pks = peaks{mm}(1:2:end);
   try
       odfVals = odf(pks, mm);
   catch % if this happens, the order of the pks has been messed up with some of the selected pks located in the second half of the odf (all the pks should be in the first half)
       pks = peaks{mm}(peaks{mm} < size(odf,1));
       odfVals = odf(pks, mm);
       fprintf('+');
   end
   %sort by odf value to arrange into fa0, fa1, etc.
   [odfVals, idx] = sort(odfVals,'descend');
   
   min_odf = min(odf(:, mm));
   
   pks = pks(idx);
   if (isempty(pks))
       pks = [1];
       odfVals = 1e-4+min_odf;
   end
   if (odfVals(1) < 1e-4+min_odf)
       odfVals(1) = 1e-4+min_odf;
   end
   for ff = 0:maxfibers-1
      if (ff+1 <= length(pks))
         % z0: see Frank Yeh, Generalized q-Space Sampling, IEEE TMI (29) 1626-35
         fal{ff+1}(mt(1), mt(2), mt(3)) = z0*(odfVals(ff+1)  - min_odf);
         indexl{ff+1}(mt(1), mt(2), mt(3)) = pks(ff+1) - 1;%NOTE: -1 to go from matlab to C++
      else
         break;
      end
   end
end

%%

%reshape the fa* and index* arrays to 1xN arrays, since v4 mat files can't
%support more than 2-d
for ff = 0:maxfibers-1
   eval(sprintf('fa%d = reshape(fal{%d}, [1 prod(dimension)]);', ff, ff+1));
   eval(sprintf('index%d = reshape(indexl{%d}, [1 prod(dimension)]);', ff, ff+1));
end %ff

%append fa* and index* variables to the fib file
save(fibfile, '-regexp','b0','index[0-9]', 'fa[0-9]', 'nshells', 'bvalue', ...
    'nmeasshell', '-mat', '-append');
if (output);fprintf(1, 'done\n');end

% make a report
report = ['A Radial Diffusion Spectrum Imaging reconstruction was applied to ' ...
    num2str(nshells) ' shell diffusion data with b-values '];
for i = 1:(length(bvalues)-1)
    report = [report num2str(bvalues(i)) ', '];
end;
if (length(bvalues) > 1) report = [report(1:(end-2)) ' and ']; end;
report = [report num2str(bvalues(end)) ' s/mm2. The shells had '];
for i = 1:(length(nmeasshell)-1)
    report = [report num2str(nmeasshell(i)) ', '];
end;
if (length(bvalues) > 1) report = [report(1:(end-2)) ' and ']; end;
report = [report num2str(nmeasshell(end)) ' diffusion samples respectively.'];
report = [report 'The inplane resolution was ' num2str(voxel_size(1)) ...
    ' x ' num2str(voxel_size(2)) ' mm with a slice thickness of ' ...
    num2str(voxel_size(3)) ' mm. Reconstruction was performed as described in ' ...
    'Baete, Yutzy and Boada. Radial q-space sampling for DSI. Magn Reson Med, ' ...
    '76(3):769-80, 2016 using Matlab-scripts available at cai2r.net -> ' ...
    'Resources -> Software Downloads.'];
report = double(report);

% see Frank Yeh, Generalized q-Space Sampling, IEEE TMI (29) 1626-35
odf = z0*odf;

%save odf* variables
if (output);fprintf(1, 'Saving ODFs...');end
odfVarNum = 0;
for start = 1:odfVarSize:size(odf,2)
  if (start-1 + odfVarSize > size(odf,2))
     eval(sprintf('odf%d=odf(:,start:end);', odfVarNum));
  else
     eval(sprintf('odf%d=odf(:,start:(start-1+odfVarSize));', odfVarNum));
  end
  odfVarNum = odfVarNum + 1;
end %start
max_odf = max(odf, [], 2)';
save(fibfile, '-regexp', 'odf[0-9]', 'max_odf', 'report', '-mat', '-append');
if (output);fprintf(1, 'done\n');end

if (output);fprintf(1, 'Compressing output file...');end
if (isunix())
    unix(['gzip -f ' fibfile ]);
else
    gzip(fibfile);
end
fibfileout = [fibfile '.gz'];
if (output);fprintf(1, 'done\n');end

%% CLEAN UP
% if the input was .src.gz, delete the temporary extracted .src
if (exist('srcfileorig', 'var'))
    delete(srcfile);
end

for i = 1:80:length(report)
    display(['  ' char(report(i:min(i+79,length(report))))]);
end;
display(['  Results saved in ' fibfileout '.']);

end %function rdsi_recon
