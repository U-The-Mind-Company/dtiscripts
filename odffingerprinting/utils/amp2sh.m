% Steven Baete
% NYU SOM CBI/CAI2R
% June 2022

% based on code by J-Donald Tournier, https://github.com/jdtournier/csd

% amp   input signals in the shape of [N,dirs]

function sh = amp2sh(amp, shtrans)

if (~isfield(shtrans,'shinv'))
    shtrans.shinv = pinv(shtrans.sh);
end;

sh = (shtrans.shinv * amp')';