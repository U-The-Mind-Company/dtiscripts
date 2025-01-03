% Steven Baete
% NYU LMC CBI
% July 2013

function v = subtractDiffDir(v1,v2)

n = size(v1,2);

v = zeros(4,n);
for i = 1:n
    if (sum(v2(1,i)) > 0)
        v(2:4,i) = v1(2:4,i)*v1(1,i) - v2(2:4,i)*v2(1,i);
        v(1,i) = sqrt(sum(v(2:4,i).^2));
        v(2:4,i) = v(2:4,i)/v(1,i);
    else
        v(:,i) = v1(:,i);
    end;
end;
