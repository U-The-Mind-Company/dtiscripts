function [Y] = sincpp(X)
%SINCPP 2nd derivative of a sinc function
%INPUTS
%  X - radians
%OUTPUTS
%  Y - 2nd derivative of a sinc evaluated at X
%
%F Boada, S Yutzy
%University of Pittsburgh
%2/2/2012

%% CALCULATIONS
Y = zeros(size(X), class(X));
%part 1: X near 0
idx = abs(X) <= 0.001;
X2 = X(idx);
Y(idx) = -1/3 + X2.*X2./10;

%part 2: everywhere else
idx = ~idx;
X2 = X(idx);

Y(idx) = 2*sin(X2)./X2./X2./X2 - 2*cos(X2)./X2./X2 - sin(X2)./X2;

end

