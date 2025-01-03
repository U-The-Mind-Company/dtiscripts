% Steven Baete
% NYU SOM CBI/CAI2R
% June 2022

% based on code by J-Donald Tournier, https://github.com/jdtournier/csd

% dirs   spherical or cartesian coordinates of the directions on the sphere

function shtrans = shtransformation(dirs, lmax)

if (size(dirs,2) == 3)
    % to spherical coordinates
    [az,el] = cart2sph(dirs(:,1),dirs(:,2),dirs(:,3)); el = pi/2 - el; % matlab has a different elevation convention
elseif (size(dirs,2) == 2)
    az = dirs(:,1);
    el = dirs(:,2);
else
    display(' shtransformation: the dirs input should be [N,2] (az,el), or [Nx3] (x,y,z)');
    return;
end;

shtrans.el = el;
shtrans.az = az;
shtrans.lmax = lmax;
sh = [];

for l = 0:2:lmax

    % lth order spherical harmonic coefficients
    s = ones(size(az,1),1);
    if (l > 0)
        s = [sqrt(2)*sin(az*(l:-1:1)) s sqrt(2)*cos(az*(1:l))];
    end

    % associated legendre polynomial 
    s2 = legendre(l, cos(el'));
    for m = 0:l
        s2(m+1,:) = s2(m+1,:).*sqrt( (2*l+1)*factorial(l-m) / ((4*pi)*factorial(l+m)));
    end
    if l
        s2 = [s2(end:-1:2,:); s2];
    end

    s = s2 .* s';

    sh = [sh s'];
end;

shtrans.sh = sh;
