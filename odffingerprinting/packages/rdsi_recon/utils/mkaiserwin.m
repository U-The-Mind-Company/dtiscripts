function w = mkaiserwin(N, beta, p, normwin)

% MKAISERWIN Symmetric Modified Kaiser-Bessel window
%    MKAISERWIN(N,BETA,P,NORMWIN) returns a symmetric N point Modified
%    Kaiser-Bessel window with shaping parameter BETA.  P defines the
%    power of the Modified Kaiser window.  When P = 1 (default), a
%    regular Modified Kaiser window is obtained.
%
%    When true (default), NORMWIN normalizes the window such that the
%    maximum amplitude is 1.
%
%    The Modified Kaiser window is also called the I1-cosh window.
%
%    See also KAISER, KAISERWIN.

%  Joe Henning - Dec 2013
        
if nargin < 2
   fprintf('??? Bad beta input to mkaiserwin ==> beta must be specified\n');
   w = [];
   return;
end

if nargin < 3
   p = 1;
end

if nargin < 4
   normwin = 1;
end

if (N == 1)
   w = 1;
   return
end

M = (N-1)/2;

w = [];
for k = 0:M
   n = k-M;
   if (k == 0)
      w(k+1) = beta/2.0/besseli(1,beta);
   else
      w(k+1) = besseli(1,beta*sqrt(1-4*n*n/(N-1)/(N-1)))/besseli(1,beta)/sqrt(1-4*n*n/(N-1)/(N-1));
   end
   w(N-k) = w(k+1);
end

% apply power
w = w.^p;

if (normwin)
   % normalize
   w = w(:)/max(w);
else
   w = w(:);
end
