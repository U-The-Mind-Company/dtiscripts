%% calculate ODFs from DWI

% Steven Baete
% NYU SOM CBI
% November 2016

function odf = calculate_ODF(dwi,F)

odf = F(1:(end/2),:)*dwi';
clear dwi;
% odf = odf(1:(end/2),:);
odf = odf';

odf(isnan(odf)) = 0;

% max_odf = max(odf,[],2);
% min_odf = min(odf,[],2);
% odf = (odf - repmat(min_odf,[1,size(odf,2)])) ...
%     ./repmat(max_odf-min_odf,[1,size(odf,2)]);
% odf = (odf - repmat(min_odf,[1,size(odf,2)]));