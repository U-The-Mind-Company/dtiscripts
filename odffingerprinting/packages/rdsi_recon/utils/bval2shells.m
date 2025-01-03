% Steven Baete
% NYU SOM CBI / CAI2R
% December 2020

function [bshells,nmeasshell,shell] = bval2shells(bvals)

[sbtable,~] = sort(bvals);
inds = find(diff(sbtable)./max(sbtable(1:(end-1)),50) > 0.1);
shellstep = [sbtable(1)/2,(sbtable(inds)+sbtable(inds+1))/2];
nmeasshell = diff([[0,inds],length(sbtable)]);
nshells = length(shellstep);
shell = zeros(size(bvals));
for i = 1:nshells
    shell(bvals >= shellstep(i)) = i;
end
for i = 1:nshells
    bshells(i) = mean(bvals(shell==i));
end;