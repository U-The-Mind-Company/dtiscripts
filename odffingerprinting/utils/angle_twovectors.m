%% calculate angle between two vectors v1 and v2

% Steven Baete
% NYU SOM CBI
% January 2017

function diffangle = angle_twovectors(v1,v2)

if (size(v1,1) == 1 & size(v2,1) > 1)
   v1 = repmat(v1,[size(v2,1),1]);
end

costheta = sum(v1.*v2,2)./(sqrt(sum(v1.^2,2)).*sqrt(sum(v2.^2,2)));
diffangle = acosd(costheta);