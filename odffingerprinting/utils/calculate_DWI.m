%% calculate DWI for voxel configurations

% Steven Baete
% NYU SOM CBI
% November 2016

%   Input:
%       vox: structure with the voxel configurations
%       q: q-matrix
%       Delta: diffusion time, can be left empty
%       delta: gradient duration, can be left empty
%   Output:
%       dwi: all dwi-values for all library elements in vox

function dwi = calculate_DWI(vox,q,Delta,delta,odf_vertices)

if ((nargin < 3) | isempty(Delta)), Delta = 1; end;
if ((nargin < 4) | isempty(delta)), delta = 0; end;
if (isfield(vox(1),'disp')) calcdispersion = true; else; calcdispersion = false; end;
if ((nargin < 5) | isempty(odf_vertices)), 
    if (calcdispersion) display(['   calculate_DWI: odf_vertices needed as'...
            ' 5th input to simulate dispersion, was not supplied!']); end;
    calcdispersion = false; end;

q = q*sqrt(1e-9);
qnorm = normalizevector(q);
nq = size(q,1);
dwi = zeros(length(vox),nq);
b = sum(q.^2,2)*(Delta-delta/3);
if (calcdispersion), odf_vertices_t = odf_vertices'; end;

if (calcdispersion)
    dispersionlist = unique(vertcat(vox(:).disp));
    if (sum(dispersionlist ~= 0) == 0) % no dispersion
        calcdispersion = false;
    end;
end;
    
if (calcdispersion)
    %% filenames
    randnum = num2str(round(rand(1)*1000));
    ampname = ['tmp' randnum '_amp.nii'];
    shname = ['tmp' randnum '_sh.nii'];
    dirname = ['tmp' randnum '_dirs_XX.nii'];
    vertname = ['tmp' randnum '_vert.dirs'];
    intpointsname = ['tmp' randnum '_intpoints_XX.nii'];
    
    %% which dispersions are needed
    dispersionlist = unique(vertcat(vox(:).disp));
    dispersionlist = dispersionlist(dispersionlist>3);
    dispersionlist(dispersionlist > 50) = 50;
    
    %% lookuptable for concentration parameters and dispersion angles
    %  according to Jelescu2016b [3]
    % psi is the angle between an individual axon and the main diffusion
    % orientation (Jelescu2017, above [5])
    % so psi is double the desired dispersion

    kappaLUT = [0.1:0.1:0.4,0.5:1:200];
    psiLUT = acos(sqrt(-1./(2*kappaLUT) ...
    +1./(sqrt(pi)*exp(-kappaLUT).*erfi(sqrt(kappaLUT)).*sqrt(kappaLUT))))*180/pi;

    %% calculate watson (simplification of bingham) distribution
    %  for each of the dispersions we need
    B = struct();
    B.d = 3; %dimension
    B.V = [0,0,1;0,1,0]';%[0,0; 0,1; 1,0]; %orthogonal direction matrix (directions of concentration)

    for k = 1:length(dispersionlist)
        kappa(k) = interp1(psiLUT,kappaLUT,dispersionlist(k)/2);
        kappa(isnan(kappa)) = 50;
        B.Z = [kappa(k),0]; % concentration parameters (second 0 for axi-symmetric)
        [B.F] = bingham_F(B.Z);

        f(k,:) = bingham_pdf(odf_vertices(:,1:(end/2))',B);
    end;
    
    % normalize these odfs
    f = f./sum(f,2)/2;
    
    f = reshape(repmat(permute(f,[3,1,2]),[nq,1,1]),[nq*length(dispersionlist),size(f,2)]);

    % save SH to mrtrix
    reamp = save_to_mrtrix_nii(f,ampname);

    % save the odf-vertices
    [azimuth,elevation,~] = cart2sph(odf_vertices(1,:),odf_vertices(2,:),odf_vertices(3,:));
    fid = fopen(vertname,'w');
    fprintf(fid,'%f %f\n',[azimuth;pi/2-elevation]);
    fclose(fid);
    
    % calculate spherical harmonics coefficients
    system(['amp2sh -quiet -force ' ...
            ' -lmax 12' ...
            ' -directions ' vertname ...
            ' ' ampname ...
            ' ' shname]);
    
    %sh = load_from_mrtrix_nii(shname,reamp);
    
    dirall = zeros(nq,2*size(odf_vertices,2)/2);
    
    for l = 1:nq 
        % negative rotate odf-vertices to q(l) -> q(l)-dirs
        ql = q(l,:);
        R(l).R = rotation_matrix_twovectors(ql'/norm(ql),[0 0 1]');
        R(l).fname = strrep(dirname,'XX',[num2str(l,'%03.0f')]);
        anglesrot = R(l).R*odf_vertices(:,1:(end/2));
        % save the rotated odf-vertices
        [azimuth,elevation,~] = cart2sph(anglesrot(1,:),anglesrot(2,:),anglesrot(3,:));
        dirall(l,:) = [azimuth,pi/2-elevation];
    end;
    dirall = repmat(dirall,[length(dispersionlist),1]);
    save_to_mrtrix_nii(dirall,dirname);
        
    % generate amplitudes in q(l)-dirs with mrtrix
    jsystem(['sh2amp -quiet -force ' ...
        ' ' shname ...
        ' ' dirname ...
        ' ' intpointsname]);
        
    intpoints = load_from_mrtrix_nii(intpointsname,reamp);
    fRF = repmat(permute(reshape(intpoints, ...
        [nq,length(dispersionlist),size(odf_vertices,2)/2]),[3,1,2]),[2,1,1]);
end;

% %% plot test
% load('odf8.mat')
% plot_odf(squeeze(fRF(:,1,:))*100,odf_vertices,odf_faces);

%%
for v = 1:length(vox)
    if (mod(v,10000) == 0)
        display(['    calculating DWI : [' num2str(v) '/' num2str(length(vox)) ']']);
    end;
    vt = vox(v);
    fibers = vt.fibers;
    microstruct = vt.microstruct;
    micro = vt.micro;
    
    St = zeros(vt.ndir,nq);
    for j=1:(vt.ndir+1)
        if (calcdispersion)
            dispersion = fibers(j,7);
        end;
        [x,y,z] = sph2cart(fibers(j,5),fibers(j,6),1);
        fdir = [x;y;z];
        
        lambda1 = fibers(j,2);
        lambda2 = fibers(j,3);
        lambda3 = fibers(j,4);
        
        d=diag([lambda1,lambda2,lambda3]);

        [x,y,z] = sph2cart(fibers(j,5),fibers(j,6),1);
        RR = rotation_matrix_twovectors([0,0,1],[x,y,z]);
        
        % The Rotated Tensor
        DD=RR*d*transpose(RR);
        
        if (calcdispersion & dispersion > 3)
            % find which fRF fits with this fiber dispersion, index k
            k = find(dispersionlist == dispersion);
            qtmp = odf_vertices_t;
        else
            qtmp = qnorm;
        end;
        
        % Calculate the Signal at the odf-vertices and b-value b
        switch microstruct
            case 'DTI'
                if (calcdispersion & dispersion > 3)
                    Stbl = fibers(j,1)*exp(-b*sum((qtmp*DD).*qtmp,2)');
                else
                    Stbl = fibers(j,1)*exp(-b.*sum((qtmp*DD).*qtmp,2));
                end;
            case 'twocomp'
                % two compartment, see Jelescu2017 eq (14)
                % (lambda_i, lambda_e, lambda_e_p, v1)
                if (calcdispersion & dispersion > 3)
                    qt = b*((qtmp*fdir).^2)';
                else
                    qt = b.*((qtmp*fdir).^2);
                end;
                Stbl = fibers(j,1) ...
                    *(micro(j,4)*exp(-micro(j,1)*qt) ...
                      + (1-micro(j,4))*exp(-micro(j,2)*qt-micro(j,3)*(b-qt)));
        end;
        
        if (calcdispersion & dispersion > 3)
            % Convolute with the fiber response function (dispersion)
            St(j,:) = sum(Stbl'.*fRF(:,:,k),1);
            % as in Anderson2005, equation 9
        else
            St(j,:) = Stbl;
        end;
    end
    % save in our output-matrix
    dwi(v,:) = sum(St,1);
end;

%% clean up temporary files
if (calcdispersion)
    system(['rm ' ampname ' ' shname ' ' strrep(dirname,'_XX','_*') ...
        ' ' vertname ' ' strrep(intpointsname,'_XX','_*')]);
end;
end
