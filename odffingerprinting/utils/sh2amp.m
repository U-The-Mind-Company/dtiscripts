% Steven Baete
% NYU SOM CBI/CAI2R
% June 2022

% based on code by J-Donald Tournier, https://github.com/jdtournier/csd

function amp = sh2amp(sh, shtrans)

amp = (shtrans.sh * sh')';