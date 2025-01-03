%% construct the rotation matrix rotating from A to B

% Steven Baete
% NYU SOM CBI
% October 2016

% see http://math.stackexchange.com/questions/180418/ ...
%   calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/  ...
%   897677#897677

function R = rotation_matrix_twovectors(A,B)

cAB = cross(A,B);

% ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
% RU = @(cAB,sAB,dAB) eye(3) + sAB + ...
%      sAB^2*(1-dAB)/(norm(cAB)^2);
%  
% R = RU(cAB,ssc(cAB),dot(A,B)); 

% faster (39.3s vs 17.7s)
sAB = [0 -cAB(3) cAB(2); cAB(3) 0 -cAB(1); -cAB(2) cAB(1) 0];
% R = eye(3) + sAB + ...
%      sAB^2*(1-dot(A,B))/(norm(cAB)^2);
R = eye(3) + sAB + ...
     sAB^2*(1-dot(A,B))/(cAB(1)^2+cAB(2)^2+cAB(3)^2);
 
if (sum(isnan(R(:))) == 9)
    R = eye(3);
end;