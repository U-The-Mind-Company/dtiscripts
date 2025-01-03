% Steven Baete
% NYU SOM CBI
% September 2017
% density compensation function, this extends h calculated above
%    the dcf is a multiplication of 3 dimensions
%       two dimensions (orthogonal to the radial lines) are each similar to the
%       2D-radial dcf, hence, ~ (qval)^2
%       the third dimension, along the radial lines, depends on the spacing
%       of the samples.

function [dcf,dcfr,dcfangle] = dcfcalc(qtable)

dimi = size(qtable,2);

qval = sqrt(sum(qtable.^2,2));
[qvalsort,ind] = sort(qval);
dcf = qval/max(qval);
if (dimi == 3)
    dcf = dcf.*dcf;
end;
dcfangle = zeros(size(dcf));
dcfr = zeros(size(dcf));
steps = qvalsort([find((qvalsort(2:end)-qvalsort(1:(end-1)))>0.1*qvalsort(1:end-1));end]);
dd = diff(steps);
dd=dd/sum([dd;dd(end)/2])*(length(dd)+1);
dcfw = 1/2*[0;dd]+1/2*[dd;dd(end)];
dcfwd = cumsum(dcfw);
stepso = [0;dcfwd(1:(end-1))]/2+[dcfwd(1:end)]/2;
for i = 1:length(steps)
    sel = (qvalsort >= 0.95*steps(i)) & (qvalsort <= 1.05*steps(i));
    dcf(ind(sel)) = dcf(ind(sel))*dcfw(i);
    if (dimi == 3)
        dcfangle(ind(sel)) = 4*pi*stepso(i).^2/sum(sel);
    else
        dcfangle(ind(sel)) = 2*pi*stepso(i)/sum(sel);
    end;
    dcfr(ind(sel)) = dcfw(i);
    if (i == 1)
        dcfr(ind(sel)) = pi*dcfw(i)^2;
        dcfangle(ind(sel)) = pi*dcfw(i)^2;
    end;
end;
