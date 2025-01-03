%% Noise estimation of dwi-data

% Steven Baete
% NYU SOM CBI
% 2017

% see Aja-Fernandez, IEEE Im Proc, 17:8, pp1383-98
%    equation 31

%   Input:
%       dwi: dwi-data
%       b_table: b-table
%   Output:
%       var_out: voxel-wise variance
%       var_gen: average variance

function [var_out,var_gen] = dwi_variance(dwi,b_table)

sdwi = size(dwi);

% nearest neighbours, mode of local sigma
v = b_table(2:4,:);
b = b_table(1,:);
% normalize dwi
minb = min(b) + (max(b)-min(b))/100;
dwi = dwi./repmat(mean((dwi(:,b < minb)),2),[1,length(b)]);
dwi(dwi(:) > 1.5) = 1.5;
% avoid /0
b(b==0) = 1;
% duplicate points
v = [v,-v];
b = [b,b];
dwi = [dwi,dwi];

% calculate nearest points for each dwi-sample
vdist = pdist2(v',v','euclidean');
shellpos = abs(repmat(b',[1,size(b,2)])-repmat(b,[size(b,2),1])) < 0.2*abs(repmat(b',[1,size(b,2)]));
vdist = vdist+(1-shellpos)*max(vdist(:))*2;
[~,ind] = sort(vdist,1);
% calculate the variance for each point
varl = zeros(sdwi);
for i = 1:sdwi(2)
    vals = dwi(:,ind(1:6,i));
    varl(:,i) = var(vals,[],2)./mean(vals,2);
end;
varl(isinf(varl)) = 0;
varl(isnan(varl)) = 0;
mv = repmat(median(varl,2),[1,size(varl,2)]);
varl = round(varl./mv*1e6).*mv/1e6;
% mode of local sigma
var_out = mode(varl,2);
var_gen = mode(varl(:));
