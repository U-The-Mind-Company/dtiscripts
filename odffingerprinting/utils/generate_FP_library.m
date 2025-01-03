%% Generate the ODF Fingerprinting library

% Steven Baete
% NYU SOM CBI
% November 2016

% clear all;
% close all;

%% generate vox for the library

time = tic;

opt.libflag = 'phantomas';

display(['  generating library for  ' opt.libflag ' ']);

bvals = [1000,2000,3000,4000];
edge = 1.25;
MAX_FIBERS = 5;
switch opt.libflag
    case 'phantomas'
        edge = 0.8;
        opt.adcsteps = 3;
        opt.microsteps = 4;
        opt.maxmicrodiff = 1.0;
        opt.maxratdiff = 1;
        opt.minangle = 0;
        [vox,libopt] = generate_directions_QA('odf8.mat',2,[],[0.7,1.5],[0.6,0.8],opt);
    case 'phantomas_3fib'
        edge = 0.8;
        opt.adcsteps = 3;
        opt.microsteps = 1;
        opt.maxmicrodiff = 0.4;
        opt.minangle = 20; % in degrees
        opt.ratiosteps = 6;
        opt.maxratdiff = 1;
        [vox,libopt] = generate_directions_QA('odf8.mat',3,[],[0.5,1.5],[0.7,0.7],opt);  
    case 'invivo'
        bvals = [200,1450,2700,4000];
        edge = 0.6;
        opt.adcsteps = 1;
        opt.microsteps = 8;
        opt.maxmicrodiff = 1;%0.4;
        opt.maxratdiff = 1;
        opt.minangle = 0;%20
        [vox,libopt] = generate_directions_QA('odf8.mat',2,[],[0.9],[0.3,1.0],opt);
    otherwise
        [vox,libopt] = generate_directions_QA('odf8.mat',2,[],[],[],opt);
end;

%% save parameters from the library

tic;
dir = NaN*ones(length(vox),MAX_FIBERS,2);
fa = NaN*ones(length(vox),MAX_FIBERS+1);
micro = NaN*ones(length(vox),MAX_FIBERS+1,size(vox(1).micro,2));
adc = NaN*ones(length(vox),MAX_FIBERS+1);
rat = NaN*ones(length(vox),MAX_FIBERS+1);
disp = NaN*ones(length(vox),MAX_FIBERS+1);
for i = 1:length(vox)
    v = vox(i);
    orig = v.directions(:,:);
    if (isempty(v.FA))
        fa(i,1) = NaN;
        adc(i,1) = v.ADC;
        dir(i,1,:) = NaN;
        rat(i,1,:) = v.fibers(1,1);
        %disp(i,1,:) = v.disp;
        micro(i,1,:) = v.micro;
    else
        fa(i,1:(v.ndir+1)) = v.FA;
        adc(i,1:(v.ndir+1)) = v.ADC;
        dir(i,1:v.ndir,:) = orig;
        rat(i,1:(v.ndir+1),:) = v.fibers(1:(v.ndir+1),1);
        disp(i,1:(v.ndir+1),:) = v.disp(1:(v.ndir+1));
        micro(i,1:(v.ndir+1),:) = v.micro(1:(v.ndir+1),:);
    end;
end;
nvox = length(vox);

display(['  save library parameters ' num2str(toc) ' s']);

%% calculate the DWI from the library
tic;

% basename = ['dsi_q_vector_rad_59_sb.txt'];
% mrtrix: dirgen -force -cartesian 90 dir90.txt
basename = ['dir90.txt'];
[q,F,odf_faces,odf_vertices] = get_qmatrix(basename,bvals*1e6,[],[],'Edge',edge,'ODFfile',libopt.anglefile);
libopt.q = q;
libopt.F = F;
libopt.bvals = bvals;
libopt.edge = edge;

dwi = calculate_DWI(vox,q,[],[],odf_vertices);

display(['  calculate dwi ' num2str(toc) ' s']);

clear vox;

%% calculate the library ODF
tic;

odf = calculate_ODF(dwi,F);
odf((odf <= 0)) = 0;

display(['  calculate odf ' num2str(toc) ' s']);
% 
% plot_odf(odf'*10,odf_vertices,odf_faces)

if (~strcmp(opt.libflag,'invivo_mc'))
    clear dwi;
end;

%% calculate the directions and QA-values traditionally 
tic;

[dirs] = find_ODF_peak(odf,odf_faces,odf_vertices);

display(['  ref values    ' num2str(toc) ' s']);

%% rotate the ODFs to the point of maximum odf-value
tic;

[odfrot,maxind,R] = rotate_ODF_to_max(odf,odf_vertices);

display(['  rotate odf    ' num2str(toc) ' s']);

clear odf;

%% rotate the directions of the library accordingly
tic;

dirrot = rotate_dir(dir,maxind,R,false);

display(['  rotate lib    ' num2str(toc) ' s']);

if (~strcmp(opt.libflag,'invivo_mc'))
    clear dir;
end;

%% save the library
tic;

% normalize the entries in the library
odfrot = normalizevector(odfrot);

if (strcmp(opt.libflag,'invivo_mc'))
    libname = ['lib.nf' num2str(libopt.max_fiber) ...
        '.' libopt.anglefile(1:(end-4)) '.n' num2str(nvox) '.' opt.libflag '.mat'];
    save(libname,'dwi','dir','libopt','fa','adc','rat','disp','micro','-v7.3');
else
    libname = ['lib.nf' num2str(libopt.max_fiber) ...
        '.' libopt.anglefile(1:(end-4)) '.n' num2str(nvox) '.' opt.libflag '.mat'];
    save(libname,'odfrot','dirrot','libopt','fa','adc','rat','disp','micro','-v7.3');
end;

display(['  saved lib as  ' libname]);
display(['  save lib      ' num2str(toc) ' s']);

display(['  total time    ' num2str(toc(time)) ' s']);
