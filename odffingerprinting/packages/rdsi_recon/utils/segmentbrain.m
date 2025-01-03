% Steven Baete
% NYU SOM CBI
% January 2016

% output:
%     maskl : coordinates of the mask points
%     mask  : whole brain mask
%     masks : masks of csf, grey and white matter

function [maskl,mask,masks] = segmentbrain(imgt,MASK_FILE,ploton,slicesel,wmslicesel,extmask)

% For spams
setenv('MKL_NUM_THREADS','1')
setenv('MKL_SERIAL','YES')
setenv('MKL_DYNAMIC','NO')

if (nargin == 1)
    writeres = false;
    basic = false;
    ploton = 0;
end;
if (nargin == 2)
    basic = false;
    ploton = 0;
end;
if (nargin > 1 && ~isempty(MASK_FILE))
    writeres = true;
else
    writeres = false;
end;
if (nargin < 4 || isempty(slicesel))
    slicesel = 1:size(imgt,3);
end;
if (nargin < 5)
    wmslicesel = 1:size(imgt,3);
end;

imgtorig = imgt;
imgt = imgt(:,:,slicesel,:);
imgtb = imgtorig(:,:,wmslicesel,:);

dimension = size(imgt);
dimorig = size(imgtorig);

%% Reorganize raw DSI-data

s = size(imgt);
sb = size(imgtb);

V = reshape(imgt,[prod(s(1:3)),s(4)])';
Vb = reshape(imgtb,[prod(sb(1:3)),sb(4)])';

% ignore all zero or obviously noise voxels
if (sum(V(1,:)==0)/size(V,2) > 0.4)
    % already masked
    masksel = V(1,:)>0;
    maskselb = Vb(1,:)>0;
else
    t = sum(abs(V),1);
    tb = sum(abs(Vb),1);
    %masksel = (t > (max(t(:))/1000));

    [~,center] = hist(t(:));
    n = mean(t(t(:) > center(2)));
    masksel = (t/n > 0.5);

    [~,center] = hist(tb(:));
    n = mean(tb(tb(:) > center(2)));
    maskselb = (tb/n > 0.5);
end;
% imagine(reshape(masksel,s(1:3)))

Vm = V(:,masksel);
Vmb = Vb(:,maskselb);

%% nnsc

param.K=3;  % learns a dictionary with 100 elements
param.numThreads=4; % number of threads

param.iter=-5;  % let us see what happens after 100 iterations.
param.lambda=0.0001*norm(single(Vm(:)));
param.verbose=false;

[H2 Wm] = nnsc(single(Vm),param);

% sort the outputs
lH2 = log(H2./repmat(H2(1,:),[size(H2,1),1]));
order = abs(mean(lH2./repmat((1:size(lH2,1))',[1,param.K]),1));
[~,ind] = sort(order,'descend');
H2=H2(:,ind);
Wm=Wm(ind,:);

Wm2 = Wm./repmat(sum(Wm,1),[param.K,1]);
Wm2(isnan(Wm2(:))) = 1;
[~,Cm] = max(Wm2,[],1);

W2 = zeros([param.K,prod(s(1:3))]);
%W2(:,masksel) = Wm;
W2(:,masksel) = Wm;
W2 = reshape(W2',[s(1:3),param.K]);

W3 = zeros([param.K,prod(s(1:3))]);
%W2(:,masksel) = Wm;
W3(:,masksel) = Wm2;
W3 = reshape(W3',[s(1:3),param.K]);

C2 = zeros([1,prod(s(1:3))]);
C2(:,masksel) = Cm;
C2 = reshape(C2',[s(1:3),1]);

%% nnsc

param.lambda=0.0001*norm(single(Vmb(:)));
param.iter=-3;  % let us see what happens after 100 iterations.

[H2b Wmb] = nnsc(single(Vmb),param);

% sort the outputs
lH2b = log(H2b./repmat(H2b(1,:),[size(H2b,1),1]));
order = abs(mean(lH2b./repmat((1:size(lH2b,1))',[1,param.K]),1));
[~,ind] = sort(order,'descend');
H2b=H2b(:,ind);
Wmb=Wmb(ind,:);

Wm2b = Wmb./repmat(sum(Wmb,1),[param.K,1]);
Wm2b(isnan(Wm2b(:))) = 1;
[~,Cmb] = max(Wm2b,[],1);

C2b = zeros([1,prod(sb(1:3))]);
C2b(:,maskselb) = Cmb;
C2b = reshape(C2b',[sb(1:3),1]);

%% general mask

% take the grey+white matter
greyw = W2(:,:,:,2)+W2(:,:,:,3);

% threshold the grey+white matter
[~,center] = hist(greyw(:),100);
n = mean(greyw(greyw(:) > center(2)));
greyw = (greyw/n > 0.5);

% fill the white matter and csf
se = strel('disk',3);
greyw = imclose(logical(greyw),se);
for i = 1:size(greyw,3)
    greywf(:,:,i) = imfill(greyw(:,:,i),'holes');
end;
mask = greywf;
if (~isempty(extmask)) % use the externally supplied mask
    mask = extmask;
end;

%% mask for CSF, grey, white matter

se = strel('disk',1); 
for i = 1:size(W2,4)
    if (i > 1)
        temp = zeros(size(mask));
        temp(C2(:).*mask(:) == i) = 1;
        % clean up csf
    else
        wmreind = find(slicesel == wmslicesel(1)):find(slicesel == wmslicesel(end));
        
        C2bt = zeros(size(mask));
        C2bt(:,:,(wmslicesel(1)-slicesel(1))+(1:length(wmslicesel))) = C2b;
        temp = zeros(size(mask));
        temp(C2bt(:).*mask(:) == i) = 1;
        if (sum(temp(:)) == 0)
            temp(C2bt(:).*mask(:) > 0) = 1;
        end;
        
        temp(:,:,[1:(wmreind(1)-1),(wmreind(end)+1):end]) = 0;
        if (sum(temp(:)) > 2000)
            se = strel('disk',1);
            temp = imerode(temp,se);
        end;        
        % let it find the blobs
        bwmap = temp == 1;
        labeledImage = bwlabeln(bwmap, 6);
        stats = regionprops(labeledImage);

        % pick the biggest (central) blob
        [msize,label] = sort([stats(:).Area],'descend');
        if (~isempty(msize))
            label = label(msize > msize(1)/10);
            % central blob
            dist = arrayfun(@(x) sqrt(sum((x.Centroid-size(temp)/2).^2)),stats(label));
            [sdist,ind] = sort(dist);
            label =  label(ind((sdist < min( sdist(1)*3,  ...
                    sqrt(sum((size(temp)/2).^2))/4) )));
            if (~isempty(label))
                temp = zeros(size(labeledImage));
                for j = 1:length(label)
                    temp = temp + (labeledImage == label(j));    
                end;
            end;
        end;
%         % only the top half of the masked brain
%         [c,ind] = max(squeeze(sum(sum(mask,1),2)));
%         temp(:,:,1:(ind-1)) = 0;
    end;
    masks{i} = temp;
end;

%% reverse the slice selection
imgt = imgtorig;
maskt = zeros(dimorig(1:3));
maskt(:,:,slicesel) = mask;mask = maskt;
for i = 1:length(masks)
    maskt = zeros(dimorig(1:3));
    maskt(:,:,slicesel) = masks{i};masks{i} = maskt;
end;
clear maskt;

%% reorganize mask data
if (length(dimorig) >= 3)
    [xt,yt,zt] = meshgrid(1:dimorig(2),1:dimorig(1),1:dimorig(3));
    maskt = reshape(mask,[1,prod(dimorig(1:3))]);
    maskl = [yt(maskt == 1);xt(maskt == 1);zt(maskt == 1)]';

    if (writeres)
        dlmwrite(MASK_FILE,mask-1,'delimiter',' '); %-1 because the file uses C++ indices
    end;
else
    [xt,yt] = meshgrid(1:dimorig(2),1:dimorig(1));
    maskt = reshape(mask,[1,prod(dimorig)]);
    maskl = [yt(maskt == 1);xt(maskt == 1)]';
end;

%% plot segmentation and masks

if ploton == 1
    temp = double(imgt(:,:,:,1));
    out = squeeze(permute(temp/max(temp(:)),[2,1,3]));
    for i = 1:param.K
        temp = zeros(dimorig(1:3));
        temp(:,:,slicesel) = W3(:,:,:,i);
        out = cat(2,out, ...
          squeeze(permute(temp/max(temp(:)),[2,1,3])));
    end;

    temp = (double(imgt(:,:,:,1)) + 1000*mask);
    out2 = squeeze(permute(temp/max(temp(:)),[2,1,3]));
    for i = 1:param.K
        temp = (double(imgt(:,:,:,1)) + 1000*masks{i});
        out2 = cat(2,out2, ...
          squeeze(permute(temp/max(temp(:)),[2,1,3])));
    end;

    imagine(cat(1,out,out2));
end;

if ploton == 2
    %slicesel=[25:58];
    temp1 = (double(imgt(:,:,:,1)) + 1000*masks{1});
    temp2 = (1000*masks{1});
    temp1 = permute(temp1,[2,1,3]);
    temp2 = permute(temp2,[2,1,3]);
    figure('Position',[0,0,1200,700],'Color','k');
    imagesc(cat(2,montageSB(temp1(:,:,wmslicesel)),montageSB(temp2(:,:,wmslicesel))));colormap('bone');
    axis off;set(gca,'Position',[0,0,0.95,0.95]);caxis([0 1200]);axis equal;
end;
