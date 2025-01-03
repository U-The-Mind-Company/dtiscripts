%% Find the peaks of ODFs using ODF Fingerprinting

% Steven Baete
% NYU SOM CBI
% November 2016

%   Input:
%       odf: odf-values
%       odf_faces, odf_vertices: odf tesselation
%       odflib: ODF FP library
%       dirslib: directions in the library
%       microlib: microstructure parameters in the library
%       MAX_FIBERS: maximum number of fibers to search for (default 5)
%       adclib: adc's in the library
%       ratlib: volume fraction of the fibers in the library
%       matchmethod: ODF Fingerprinting option
%           0: without fiber complexity penalty
%           7: with fiber complexity penalty (default), needs var-input
%       var: variance of DWI data
%   Output:
%       dirs: identified fiber directions
%       fpmicro: identified qa-values
%       fpadc: identified adc-values
%       fprat: identified fiber volume fractions
%       fpdisp: identified fiber dispersions
function [dirs,fpmicro,fpadc,fprat,fpdisp,fprmse,fpind] = find_ODF_peak_FP(odf,odf_vertices,odf_faces,odflib,...
    dirslib,microlib,MAX_FIBERS,adclib,ratlib,matchmethod,var,displib)

if (nargin < 7 || isempty(MAX_FIBERS))
    MAX_FIBERS = 5;
end
if (nargin < 10 || isempty(matchmethod))
    matchmethod = 7;
end
if (nargin < 12 || isempty(displib))
    displib = zeros(size(adclib))*NaN;
end
if (matchmethod == 7 && (nargin < 11 || isempty(var)))
    display([' find_ODF_peak_FP: for the penalty correcting for the noise, ' ...
        ' a variance measure is needed ']);
    error([' find_ODF_peak_FP: for the penalty correcting for the noise, ' ...
        ' a variance measure is needed ']);
end

%matlab version, i.e. can we use implicit expansion?
%if newer than > R2016b
tf = ~(isMATLABReleaseOlderThan("R2016b"));

[odfrot,maxind,R] = rotate_ODF_to_max(odf,odf_vertices);

MAX_FIBERIND = min(size(microlib,2),MAX_FIBERS+1);
fpdirrot = zeros(size(odfrot,1),MAX_FIBERS,2);
fpmicro = zeros(size(odfrot,1),MAX_FIBERS+1,size(microlib,3));
fprat = zeros(size(odfrot,1),MAX_FIBERS+1);

odfrot = normalizevector(odfrot - min(odfrot,[],2));  
odflib = normalizevector(odflib - min(odflib,[],2)); 

n = size(odfrot,2);
switch(matchmethod)
    case 0 % dot-product matching method
        [~,fpind] = max(odflib(:,:)*odfrot');
    case 7 % dot-product with penalty for nr of fibers fitted 
           % (Akaike Information Criterion) and noise-correction
        pars = 1 + 4*sum(~isnan(microlib(:,:,1)),2) + sum(~isnan(displib(:,:,1)),2);
        batchsize = 30000;
        fpind = [];maxm = [];
        for k = 1:batchsize:size(odflib,1)
            kind = k-1+(1:batchsize);kind = kind(kind <= size(odflib,1));
            if (tf)
                [maxt,fpindt] = max(log(odflib(kind,:)*odfrot') ...
                        - (pars(kind) / (4*n) ) ...
                        .* max((sqrt(var')-0.003)/(0.015-0.003),0)); 
            else
                [maxt,fpindt] = max(log(odflib(kind,:)*odfrot') ...
                        -repmat(pars(kind) / (4*n),[1,size(odfrot,1)]) ...
                        .* max((repmat(sqrt(var'), ...
                        [length(kind),1])-0.003)/(0.015-0.003),0)); 
            end
            [maxm,ind2] = max([maxm;maxt],[],1);
            fpindt = [fpind;fpindt+k-1];
            fpind = fpindt(sub2ind(size(fpindt),ind2,1:size(fpindt,2)));
        end
end

fprmse = sqrt(sum((odfrot - odflib(fpind,:)).^2,2)/n);

fpdirrot(1:size(odfrot,1),1:MAX_FIBERS,:) = dirslib(fpind,1:MAX_FIBERS,:);
fpmicro(1:size(odfrot,1),1:MAX_FIBERIND,:) = microlib(fpind,1:MAX_FIBERIND,:);
if (nargin >= 8 && ~isempty(ratlib))
    fprat = zeros(size(odfrot,1),MAX_FIBERS+1);
    fprat(1:size(odfrot,1),1:(MAX_FIBERS+1)) = ratlib(fpind,1:(MAX_FIBERS+1));
end
if (nargin >= 12 && ~isempty(displib))
    fpdisp = zeros(size(odfrot,1),MAX_FIBERS+1);
    fpdisp(1:size(odfrot,1),1:(MAX_FIBERS+1)) = displib(fpind,1:(MAX_FIBERS+1));
else
    fpdisp = [];
end
if (nargin >= 7 && ~isempty(adclib))
    fpadc = zeros(size(odfrot,1),MAX_FIBERS+1);
    fpadc(1:size(odfrot,1),1:(MAX_FIBERS+1)) = adclib(fpind,1:(MAX_FIBERS+1));
else
    fpadc = [];
end

dirs = rotate_dir(fpdirrot,maxind,R,true);

end