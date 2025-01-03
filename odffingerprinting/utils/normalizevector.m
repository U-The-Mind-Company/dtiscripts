%% Normalizes vectors row by row

% Steven Baete
% NYU SOM CBI
% Oktober 2016

function [x,n] = normalizevector(x)

for i = 1:size(x,1)
    nt = norm(x(i,:));
    if (nt~=0)
        x(i,:) = x(i,:)/nt;
    end;
    n(i) = nt;
end;